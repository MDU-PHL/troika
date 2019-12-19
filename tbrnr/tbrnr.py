"""package_name 
.. moduleauthor:: Kristy Horan <kristyhoran15@gmail.com>


"""


import logging
import argparse
import configargparse
import pathlib
import sys
import os
from package_name.Class_file_1 import SomeClass

def run_pipeline(args):
    '''
    Run the pipeline for the first time
    '''
    R = SomeClass(args)
    return(R.run_pipeline())


def main():
    # setup the parser
  
    parser = configargparse.ArgumentParser(description='package_name - Some description',formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    # subparser for running the pipeline
    # if using sub commands use add_subparser otherwise 
    subparsers = parser.add_subparsers(help="Task to perform") # remove if not using sub commands
    parser_sub_run = subparsers.add_parser('run', help='Initial run of Bohra', formatter_class=configargparse.ArgumentDefaultsHelpFormatter,default_config_files=[f"{pathlib.Path.cwd().absolute() / 'bohra.conf'}"])
    
    # use add_argument
    parser.add_argument('--longform', '-l', help='Help description', default = 'default_value')

    
    parser.set_defaults(func=run_pipeline)
    args = parser.parse_args()
    
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        
        args.func(args)
	
if __name__ == '__main__':
    main()

