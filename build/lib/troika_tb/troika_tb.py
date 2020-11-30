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
from troika_tb.RunTroika import Troika
from troika_tb.versions import db_version, troika_version


def run_pipeline(args):
    '''
    Run the pipeline
    '''
    T = Troika(args)
    return(T.run_pipeline())

def main():
    # setup the parser
  
    parser = configargparse.ArgumentParser(description='Troika - a pipeline for phylogenentic analysis, detection and reporting of genomic AST in Mtb',formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + troika_version)
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
        '--profiler_singularity_path',
        '-ps',
        help='URL for TB-profiler singularity container.',
        default = 'docker://mduphl/mtbtools'
    )
    parser.add_argument(
        '--snippy_singularity_path',
        '-ss',
        help='URL for Snippy singularity container.',
        default = 'docker://mduphl/snippy:v4.4.3'

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
    parser.add_argument(
        '--profiler_threads',
        '-t',
        help='Number of threads to run TB-profiler', 
        default=1
    )
    parser.add_argument(
        '--kraken_threads',
        '-kt',
        help = 'Number of threads for kraken',
        default = 4
    )
    parser.add_argument(
        '--kraken_db', 
        '-k', 
        env_var='KRAKEN2_DEFAULT_DB', 
        help='Path to DB for use with kraken2, if no DB present speciation will not be performed.'
    )
    parser.add_argument(
        '--snippy_threads',
        '-st',
        help = 'Number of threads for snippy',
        default = 8
    )
    parser.add_argument(
        '--mode',
        '-m',
        help = "If running for MDU service use 'mdu', else use 'normal'",
        default = 'normal',
        choices = ['mdu', 'normal']
    )
    parser.add_argument(
        '--positive_control',
        '-pc',
        help = 'Path to positive control - REQUIRED if running for MDU service',
        default= '')

    parser.add_argument(
        '--db_version',
        help=f'The version of database being used.', 
        default=db_version
    )
    parser.add_argument(
        '--min_cov',
        '-mc',
        default=40,
        help = f"Minimum coverage for quality checks, isolates with coverage below this threshold will not be used in the analysis."
    )
    parser.add_argument(
        '--min_aln',
        '-ma    ',
        default=80,
        help = f"Minimum alignment for phylogenetic analysis, alignments lower than this threshold will not be included in the calculation of core-genome."
    )
    parser.set_defaults(func=run_pipeline)
    args = parser.parse_args()
    
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        
        args.func(args)
	
if __name__ == '__main__':
    main()

