import pathlib
import os, getpass, shutil, re, psutil
import pandas
import jinja2
import sh
import logging
import filecmp
import datetime
import numpy
import itertools
import subprocess
import json
from Bio import SeqIO, Phylo
from packaging import version
from troika.versions import db_version



class Troika(object):
    '''
    A class for Troika
    '''
    
    def __init__(self, args):
        # from logger import logger
        # get date and time
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        # create file handler which logs even debug messages
        fh = logging.FileHandler('troika.log')
        fh.setLevel(logging.INFO)
        # create console handler with a higher log level
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        # create formatter and add it to the handlers
        formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
        self.logger.addHandler(fh)
        self.logger.info(f"Initialising Troika")
        self.now = datetime.datetime.today().strftime("%d_%m_%y_%H")
        self.date = datetime.datetime.today().strftime("%d_%m_%y")
        self.input_file = self.check_input_file(args.input_file)
        self.workdir = self.check_input_file(args.workdir)
        self.resources = self.check_input_file(args.resources)
        self.run_singulairty = args.Singularity
        self.singularity_path = args.singularity_path
        self.resistance_only = args.resistance_only
        self.detect_species = args.detect_species
        self.db_version = db_version
        self.jobs = args.jobs  
        self.profiler_threads = int(args.profiler_threads)
        self.kraken_threads = int(args.kraken_threads)
        self.kraken_db = args.kraken_db
        self.snippy_threads = int(args.snippy_threads)
        self.profiler_version = db_version
        self.isolates = ''
        # self.output = args.output
        self.set_snakemake_jobs()
        self.min_cov = args.min_cov
        self.min_aln = args.min_aln
        self.mode = 'mdu' if args.mode else 'normal'
        if self.mode == 'mdu':
            self.check_positive_control(args.positive_control)
        else:
            self.positive_control = {}


    def check_positive_control(self, positive_control):
        if positive_control == '' or not pathlib.Path(positive_control).exists():
            self.logger.warning(f"You are running Troika for MDU service. You must include a path positive control.")
            raise SystemExit
        else:
            p = pathlib.Path(positive_control)
            reads = sorted(p.glob(f"*.f*q.gz"))
            if len(reads) == 2:
                self.positive_control = {0:"2999-99887", 1:reads[0], 2:reads[1]}
                # return pos_string
            else:
                self.logger.warning(f"The path to positive control has not been found. Please check your settings and try again.")
                raise SystemExit
            

    def set_snakemake_jobs(self):
        '''
        set the number of jobs to run in parallel based on the number of cpus from args
        '''
        self.logger.info(f"Determining number of jobs to run in parallel.")
        threads = max([self.kraken_threads, self.profiler_threads, self.snippy_threads])
        # if int(self.jobs * threads) > int(psutil.cpu_count()):
        #     self.jobs = int(psutil.cpu_count()) - 1
    
        self.logger.info(f"Troika will run : {self.jobs} in parallel... keeping it nice.")
    
    def three_cols(self, tab):
        '''
        Ensure that there are 3 columns, isolate, R1 and R2
        returns True if 3 columns False otherwise
        
        '''
        self.logger.info(f"Checking that input file is the correct structure.")
        if tab.shape[1] == 3:
            return True
        else:
            return False

    def all_data_filled(self, tab):
        '''
        Ensure that all fields contain data - no NA's
        returns True if there are no nan, False otherwise
        '''
        self.logger.info("Checking that there is no empty fields in the input file.")
        return tab.isnull().sum().sum() == 0
    
    def link_reads(self, read_source, isolate_id, r_pair):
        '''
        check if read source exists if so check if target exists - if not create isolate dir and link. If already exists report that a dupilcation may have occured in input and procedd

        '''
        # check that job directory exists
        self.logger.info(f"Checking that reads are present.")
        # check that READS exists
        R = self.workdir / isolate_id
        if not R.exists():
            R.mkdir()
        
        if f"{read_source}"[0] != '/':
            read_source = self.workdir / read_source
        
        if read_source.exists():
            read_target = R / f"{r_pair}"
            if not read_target.exists():
                read_target.symlink_to(read_source)
        else:
            self.logger.warning(f"{read_source} does not seem to a valid path. Please check your input and try again.")
            raise SystemExit()
    
    def path_exists(self,path, v = True):
        '''
        ensure files are present, if so continues if not quits with FileNotFoundError
        input:
            :path: patht to files for pipeline
            :v: if v == True print message, else just check
        output:
            returns True (or fails with FileNotFoundError)
        '''
        
        if not path.exists():
            self.logger.warning(f"The {path.name} does not exist.")
            raise FileNotFoundError(f"{path.name}")
        else:
            if v == True:
                self.logger.info(f"Found {path.name}.")

            return True

    def check_reads_exists(self, tab):
        '''
        check that the correct paths have been given
        if reads not present path_exists will cause a FileNotFound error and warn user
        :input
            :tab: dataframe of the input file
        '''
        self.logger.info(f"Checking that all the read files exist.")
        for i in tab.itertuples():
            
            if not '#' in i[1]:
                r1 = i[2]
                r2 = i[3]
                self.path_exists(pathlib.Path(r1), v = False)
                self.link_reads(pathlib.Path(r1), isolate_id=f"{i[1].strip()}", r_pair='R1.fq.gz')
                self.path_exists(pathlib.Path(r2), v = False)
                self.link_reads(pathlib.Path(r2), isolate_id=f"{i[1].strip()}", r_pair='R2.fq.gz')
        return True


    def get_samples(self):
        '''
        get samples list after ensuring that structure of input is correct and reads can be found
        '''

        tab = pandas.read_csv(self.input_file, sep = '\t', header = None)
        if not self.three_cols(tab):
            logging.warning(f"{self.input_file} does not appear to be in the correct configuration")
            raise TypeError(f"{self.input_file} has incorrect number of columns")
        
        if not self.all_data_filled(tab):
            self.logger.warning('warning',f"{self.input_file} appears to be missing some inforamtion.")
            raise TypeError(f"{self.input_file} appears to be missing some inforamtion.")
        
        if tab.shape[0] < 4:
            self.resistance_only = True
        
        if self.mode == 'mdu':
            dx = pandas.DataFrame(self.positive_control, index = [0])
            tab = tab.append(dx, sort = True)

        self.check_reads_exists(tab = tab)
        
        self.isolates = ' '.join(list(tab.iloc[:,0]))
    
    def check_input_file(self, path):
        '''
        Check that all files and directories exist
        '''
        self.logger.info(f"Searching for : {path}")
        if pathlib.Path(path).exists():
            self.logger.info(f"Success : {path} found.")
            return pathlib.Path(path)
        else:
            self.logger.warning(f"Failed : {path} was not found, please check your input and try again.")
            raise SystemExit
    
    def tbprofiler(self):
        '''
        check the installation of tb-profiler
        '''
        self.logger.info(f"Checking that tb-profiler is installed and recording version.")
        version_pat = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+))?\b')
        try:
            tp = subprocess.run(['tb-profiler', 'version'], capture_output = True, encoding='utf-8')
            tp = tp.stdout
            # self.logger.info(f"{tp}")
            self.profiler_version = version_pat.search(tp).group(0)
            self.logger.info(f"TB-Profiler version {self.profiler_version} found. Good job!")
    
        except FileNotFoundError:
            self.logger.warning(f"TB-Profiler is not installed, please read installation/usage instructions and try again.")
            raise SystemExit

    def check_installation(self):
        '''
        If not using singularity then check for tbprofiler installation
        '''
        
        if not self.run_singulairty:
            self.tbprofiler()
        else:
            self.logger.info(f"You are using a singularity container from {self.singularity_path} today.")

    def generate_smk_config(self):
        '''
        Generate the config file
        '''
        config = {
                'template_path': f"{self.resources / 'templates'}",
                'script_path':f"{self.resources / 'utils'}", 
                'samples':self.isolates,
                'singularity_path_profiler' : self.singularity_path, 
                'profiler_threads': self.profiler_threads, 
                # 'final_output' : final_output, 
                'db_version': self.db_version,
                'run_species': self.detect_species,
                'kraken_db': self.kraken_db,
                'kraken_threads':self.kraken_threads, 
                'snippy_threads': self.snippy_threads,
                'amr_only': self.resistance_only, 
                'reference': f"{self.resources / 'reference'/ 'tbdb.fasta'}", 
                'index': f"{self.resources / 'reference'/ 'tbdb.fasta.fai'}", 
                'mask': f"{self.resources / 'reference'/ 'mask.bed'}",
                'mode': self.mode,
                'positive_control': self.positive_control,
                'min_cov':self.min_cov,
                'min_aln':self.min_aln
        }

        config_source = self.resources / "templates" / "config.yaml"
        self.logger.info(f"Writing config file")
        config_template = jinja2.Template(config_source.read_text())
        config_target = self.workdir / "config.yaml"
        config_target.write_text(
            config_template.render(
                config
            )
        )
  
    def construct_command(self):
        '''
        once all checks have passed then construct the command
        '''
        # TODO add support for sbatch and qsub
        singularity = f"--use-singularity --singularity-args '--bind /home'" if self.run_singulairty else ''
        cmd = f"snakemake -s {self.resources / 'utils' / 'troika.smk'} -j {self.jobs} -d {self.workdir} {singularity}"
        return cmd
        

    def run_snakemake(self):
        '''
        run the snakemake command
        '''
        cmd = self.construct_command()
        self.logger.info(f"Now running AMR detection using the Troika pipeline with command {cmd}. Please be patient this may take some time.")
        wkf = subprocess.run(cmd, shell = True, capture_output = True)
        if wkf.returncode == 0:
            return True
        else:
            return False


    def run_pipeline(self):
        '''
        run the pipeline
        '''
        
        # run checks
        self.logger.info(f"Checking TB-profiler installation.")
        self.check_installation()
        # check input file and get samples
        self.logger.info(f"Checking input file structure is correct and obtaining a list of samples.")
        self.get_samples()
        self.logger.info(f"Generating a config file for todays Troika job.")
        self.generate_smk_config()
        self.logger.info(f"All checks and preparations are completed.")
        if self.run_snakemake():
            self.logger.info(f"Checking that all outputs have been generated.")
            if self.mode == 'normal' and pathlib.Path(f"troika.tab").exists():
                self.logger.info(f"SUCCESS - all outputs have been generated. Have a nice day!")
            elif self.mode == 'mdu' and pathlib.Path(f"MMS155_{self.date}.csv").exists():
                self.logger.info(f"SUCCESS - all outputs have been generated. Have a nice day!")
            else:
                self.logger.warning(f"Something has gone wrong, output files have not been correcty created. Please check the logs and try again.")
    