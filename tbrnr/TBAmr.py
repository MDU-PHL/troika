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
from logger import logger
import versions



class TBAmr(object):
    '''
    A class for TBAmr
    '''
    
    def __init__(self, args):
        # get date and time
        logger.info(f"Initialising TBAmr")
        self.now = datetime.datetime.today().strftime("%d_%m_%y_%H")
        self.day = datetime.datetime.today().strftime("%d_%m_%y")
        self.input_file = self.check_input_file(args.input_file)
        self.workdir = self.check_input_file(args.workdir)
        self.resources = self.check_input_file(args.resources)
        self.run_singulairty = args.Singularity
        self.singularity_path = f"shub://phgenomics-singularity/tbprofiler" if args.singularity_path == '' else self.check_input_file(args.singularity_path)
        self.db_version = versions.db_version
        self.jobs = args.jobs  
        self.threads = args.threads
        self.profiler_version = versions.db_version
        self.isolates = ''
        self.output = args.output
        self.set_snakemake_jobs()

    def set_snakemake_jobs(self):
        '''
        set the number of jobs to run in parallel based on the number of cpus from args
        '''
        logger.info(f"Determining number of jobs to run in parallel.")

        if int(self.jobs * self.threads) > int(psutil.cpu_count()):
            self.jobs = int(psutil.cpu_count()) - 1
    
        logger.info(f"TBrnr will run : {self.jobs} in parallel... keeping it nice.")
    
    def three_cols(self, tab):
        '''
        Ensure that there are 3 columns, isolate, R1 and R2
        returns True if 3 columns False otherwise
        
        '''
        logger.info(f"Checking that input file is the correct structure.")
        if tab.shape[1] == 3:
            return True
        else:
            return False

    def all_data_filled(self, tab):
        '''
        Ensure that all fields contain data - no NA's
        returns True if there are no nan, False otherwise
        '''
        logger.info("Checking that there is no empty fields in the input file.")
        return tab.isnull().sum().sum() == 0
    
    def link_reads(self, read_source, isolate_id, r_pair):
        '''
        check if read source exists if so check if target exists - if not create isolate dir and link. If already exists report that a dupilcation may have occured in input and procedd

        '''
        # check that job directory exists
        logger.info(f"Checking that reads are present.")
        # check that READS exists
        R = self.workdir / 'READS'
        if not R.exists():
            R.mkdir()
        
        if f"{read_source}"[0] != '/':
            read_source = self.workdir / read_source
        
        if read_source.exists():
            I = R / f"{isolate_id}" # the directory where reads will be stored for the isolate
            if not I.exists():
                I.mkdir()
            read_target = I / f"{r_pair}"
            if not read_target.exists():
                read_target.symlink_to(read_source)
        else:
            logger.warning(f"{read_source} does not seem to a valid path. Please check your input and try again.")
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
            logger.warning(f"The {path.name} does not exist.")
            raise FileNotFoundError(f"{path.name}")
        else:
            if v == True:
                logger.info(f"Found {path.name}.")

            return True

    def check_reads_exists(self, tab):
        '''
        check that the correct paths have been given
        if reads not present path_exists will cause a FileNotFound error and warn user
        :input
            :tab: dataframe of the input file
        '''
        logger.info(f"Checking that all the read files exist.")
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
            logger.warning('warning',f"{self.input_file} appears to be missing some inforamtion.")
            raise TypeError(f"{self.input_file} appears to be missing some inforamtion.")
        
        self.check_reads_exists(tab = tab)
        
        self.isolates = ' '.join(list(tab.iloc[:,0]))
    
    def check_input_file(self, path):
        '''
        Check that all files and directories exist
        '''
        logger.info(f"Searching for : {path}")
        if pathlib.Path(path).exists():
            logger.info(f"Success : {path} found.")
            return pathlib.Path(path)
        else:
            logger.warning(f"Failed : {path} was not found, please check your input and try again.")
            raise SystemExit
    
    def tbprofiler(self):
        '''
        check the installation of tb-profiler
        '''
        logger.info(f"Checking that tb-profiler is installed and recording version.")
        version_pat = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+))?\b')
        try:
            tp = subprocess.run(['tb-profiler', 'version'], capture_output = True, encoding='utf-8')
            tp = tp.stdout
            # logger.info(f"{tp}")
            self.profiler_version = version_pat.search(tp).group(0)
            logger.info(f"TB-Profiler version {self.profiler_version} found. Good job!")
    
        except FileNotFoundError:
            logger.warning(f"TB-Profiler is not installed, please read installation/usage instructions and try again.")
            raise SystemExit

    def check_installation(self):
        '''
        If not using singularity then check for tbprofiler installation
        '''
        
        if not self.run_singulairty:
            self.tbprofiler()
        else:
            logger.info(f"You are using a singularity container from {self.singularity_path} today.")

    def generate_smk_config(self):
        '''
        Generate the config file
        '''
        final_output = f"{self.output}.csv,{self.output}.json"
        config_source = self.resources / "templates" / "config.yaml"
        logger.info(f"Writing config file")
        config_template = jinja2.Template(config_source.read_text())
        config_target = self.workdir / "config.yaml"
        config_target.write_text(
            config_template.render(
                script_path=f"{self.resources / 'utils'}", samples=self.isolates,singularity_path = self.singularity_path, threads = self.threads, final_output = final_output, db_version = self.db_version
            )
        )
    
    def construct_command(self):
        '''
        once all checks have passed then construct the command
        '''
        # TODO add support for sbatch and qsub
        cmd = f"snakemake -s {self.resources / 'templates' / 'tbamr.smk'} -j {self.jobs} -d {self.workdir} 2>&1 | tee -a {self.workdir / 'tbrnr.log'}"
        return cmd
        # snakemake -s "/opt/conda/lib/python3.7/site-packages/abritamr/templates/Snakefile.smk" -j 16 -d /home/khhor/MDU/JOBS/ad_hoc/abritamr_in_parallel/20191219  2>&1 | tee -a /home/khhor/MDU/JOBS/ad_hoc/abritamr_in_parallel/20191219/job.log.

    def run_snakemake(self):
        '''
        run the snakemake command
        '''
        cmd = self.construct_command()
        logger.info(f"Now running AMR detection using the TBrnr pipeline with command {cmd}. Please be patient this may take some time.")
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
        logger.info(f"Checking TB-profiler installation.")
        self.check_installation()
        # check input file and get samples
        logger.info(f"Checking input file structure is correct and obtaining a list of samples.")
        self.get_samples()
        logger.info(f"Generating a config file for todays TBrnr job.")
        self.generate_smk_config()
        logger.info(f"All checks and preparations are completed.")
        if self.run_snakemake():
            logger.info(f"Checking that all outputs have been generated.")
            if pathlib.Path(f"{self.output}.csv").exists() and pathlib.Path(f"{self.output}.json").exists():
                logger.info(f"SUCCESS - all outputs have been generated. Have a nice day!")
            else:
                logger.warning(f"Something has gone wrong, output files have not been correcty created. Please check the logs and try again.")
    