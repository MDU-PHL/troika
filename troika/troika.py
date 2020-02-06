"""troika
.. moduleauthor:: Kristy Horan <kristyhoran15@gmail.com>


"""


import logging
import argparse
import configargparse
import pathlib
import sys
import os
import datetime
import RunTroika
import versions

def run_pipeline(args):
    '''
    Run the pipeline
    '''
    T = RunTroika.Troika(args)
    return(T.run_pipeline())

def main():
    # setup the parser
  
    parser = configargparse.ArgumentParser(description='Troika - a pipeline for phylogenentic analysis, detection and reporting of genomic AST in Mtb',formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    # use add_argument
    parser.add_argument('--input_file',
        '-i', 
        help='Input file tab-delimited file 3 columns isolate_id path_to_r1 path_to_r2',
        default = ''
    )
    parser.add_argument('--detect_species',
        '-d',
        action = 'store_true',
        help = 'Set if you would like to detect species - note if not set troika may include non-tuberculosis species in the analysis.'
    )
    parser.add_argument('--resistance_only',
        action = 'store_true',
        help = 'If detection of resistance mutations only is needed. Phylogeny will not be performed.'
    )
    parser.add_argument(
        "--Singularity",
        "-S",
        action="store_true",
        help="If singularity is to be used for troika."
    )
    parser.add_argument(
        '--singularity_path',
        '-s',
        help='URL for TB-profiler singularity container.',
        default = 'docker://mduphl/tbprofiler'
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
        "--jobs", "-j", default=8, help="Number of jobs to run in parallel."
    )
    parser.add_argument('--profiler_threads',
        '-t',
        help='Number of threads to run TB-profiler', 
        default=1
    )
    parser.add_argument('--kraken_threads',
        '-kt',
        help = 'Number of threads for kraken',
        default = 4
    )
    parser.add_argument('--kraken_db', 
        '-k', 
        env_var='KRAKEN2_DEFAULT_DB', 
        help='Path to DB for use with kraken2, if no DB present speciation will not be performed.'
    )
    parser.add_argument('--snippy_threads',
        '-st',
        help = 'Number of threads for snippy'
    )
    parser.add_argument('--output',
        '-o',
        help = "Name of output files",
        default=f"tbrnr_{datetime.datetime.today().strftime('%d_%m_%y')}"
    )
    parser.add_argument('--db_version',
        help=f'The version of database being used.', 
        default=versions.db_version
    )

    parser.set_defaults(func=run_pipeline)
    args = parser.parse_args()
    
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        
        args.func(args)
	
if __name__ == '__main__':
    main()

