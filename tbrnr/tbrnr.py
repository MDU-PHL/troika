"""package_name 
.. moduleauthor:: Kristy Horan <kristyhoran15@gmail.com>


"""


import logging
import argparse
import configargparse
import pathlib
import sys
import os
from tbrnr.TBAmr import TBAmr

def run_pipeline(args):
    '''
    Run the pipeline
    '''
    T = TBAmr(args)
    return(T.run_pipeline())


def main():
    # setup the parser
  
    parser = configargparse.ArgumentParser(description='TBrnr - a pipeline for detection and reporting of genomic AST in Mtb',formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    # use add_argument
    parser.add_argument('--input_file',
        '-i', 
        help='Input file tab-delimited file 3 columns isolate_id path_to_r1 path_to_r2',
        default = '')
    parser.add_argument(
        "--Singularity",
        "-S",
        action="store_true",
        help="If using singularity container for AMRfinderplus"
    )
    parser.add_argument(
        "--singularity_path",
        "-s",
        default="",
        help="Path to the singularity container for TB-profiler, if empty will defualt to shub://phgenomics-singularity/tbprofiler"
    )
    parser.add_argument(
        "--workdir",
        "-w",
        default=f"{pathlib.Path.cwd().absolute()}",
        help="Working directory, default is current directory",
    )
    parser.add_argument(
        "--resources",
        "-r",
        default=f"{pathlib.Path(__file__).parent }",
        help="Directory where templates are stored",
    )
    parser.add_argument(
        "--jobs", "-j", default=16, help="Number of jobs to run in parallel."
    )

    parser.add_argument('--cpus',
        '-c',
        help='Number of CPU cores to run, will define how many rules are run at a time', 
        default=36
    )
    parser.set_defaults(func=run_pipeline)
    args = parser.parse_args()
    
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        
        args.func(args)
	
if __name__ == '__main__':
    main()

