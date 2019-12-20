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
from package_name.logger import logger



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
        
        self.jobs = args.jobs  
        self.cpus = args.cpus 

        self.set_snakemake_jobs()

    def set_snakemake_jobs(self):
        '''
        set the number of jobs to run in parallel based on the number of cpus from args
        '''
        logger.info(f"Determining number of jobs to run in parallel.")
        if int(self.cpus) < int(psutil.cpu_count()):
            self.jobs =  self.cpus
        else:
            self.jobs = 1
        logger.info(f"TBrnr will run : {self.job} in parallel... keeping it nice.")
    
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
        

    def check_singularity_path(self):
        '''
        if using singularity check the path
        '''
        if self.singularity_path != f"shub://phgenomics-singularity/tbprofiler":


    def check_installation(self):
        '''
        If not using singularity then check for tbprofiler installation
        '''
        pass

    def generate_smk_config(self):
        '''
        Generate the config file
        '''
        pass
    
    def construct_command(self):
        '''
        once all checks have passed then construct the command
        '''
        pass

    def run_snakemake(self):
        '''
        run the snakemake command
        '''
        pass

    def run_pipeline(self):
        '''
        run the pipeline
        '''
        pass
    
    