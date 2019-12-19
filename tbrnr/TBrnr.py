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



class SomeClass(object):
    '''
    A class for some_package
    '''
    
    def __init__(self, args):
        # get date and time
        self.now = datetime.datetime.today().strftime("%d_%m_%y_%H")
        self.day = datetime.datetime.today().strftime("%d_%m_%y")
        

    def some_function(self):
        '''
        some description
        '''
        pass
    
    