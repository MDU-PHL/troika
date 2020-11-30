import pandas, pathlib, subprocess

configfile:'config.yaml'

SAMPLE = config['samples'].split()

SINGULARITY_PATH_PROFILER = config['singularity_path_profiler'] # path to container for running tb-profiler
SINGULARITY_PATH_SNIPPY = config['singularity_path_profiler'] # path to container for running tb-profiler
PROFILER_THREADS = config['profiler_threads']
SCRIPT_PATH = config['script_path']
 # a two (or three item) item list - jobid.csv, jobid.json, mtb.tab
MODE = config['mode']
POSITIVE_CONTROL = config['positive_control']
DB_VERSION = config['db_version']
RUN_SPECIES = config['run_species']
KRAKEN_DB = config['kraken_db']
# PREFILL_PATH = config['prefill_path']
KRAKEN_THREADS = config['kraken_threads']
SNIPPY_THREADS = config['snippy_threads']
REFERENCE = config['reference']
IDX = config['index']
MASK = config['mask']
AMR_ONLY = config['amr_only']
WORKDIR = f"'{pathlib.Path.cwd().absolute()}'"
TEMPLATE_PATH = config['template_path']
MIN_COV = config['min_cov']
MIN_ALN = config['min_aln']

rule all:
    input:
        # REPORT_INPUT,
        expand("{sample}/snippy.toml", sample = SAMPLE),
        'kraken.toml', 'iqtree.toml', 'report.toml', 'seqdata.toml', 'distances.toml','snippy_core.toml','resistance.toml', expand("{sample}/snpit.toml", sample = SAMPLE)

rule run_kraken:
    input:
        r1='{sample}/R1.fq.gz',
        r2='{sample}/R2.fq.gz'
    output:
        "{sample}/kraken.toml"
    params:
        kraken_db = KRAKEN_DB,
        script_path = SCRIPT_PATH,
        run_kraken = RUN_SPECIES
    script:
        "kraken.py"

rule combine_kraken:
    input:
        expand("{sample}/kraken.toml", sample = SAMPLE)
    output:
        "kraken.toml"
    params:
        script_path = SCRIPT_PATH
    script:
        "combine_kraken.py"

rule estimate_coverage:
	input:
		r1="{sample}/R1.fq.gz",
		r2="{sample}/R2.fq.gz"
	output:
		"{sample}/mash.toml"
	params:
		script_path = SCRIPT_PATH
	
	script:
		"mash.py"

rule seqdata:
	input:
		r1 = '{sample}/R1.fq.gz',
		r2 = '{sample}/R2.fq.gz',
		mash = '{sample}/mash.toml'
	output:
		"{sample}/seqdata.toml"
	params:
		script_path=SCRIPT_PATH,
		mincov = MIN_COV
	
	script:
		"seqdata.py"

rule combine_seqdata:
	input:
		expand("{sample}/seqdata.toml", sample = SAMPLE)
	output:
		"seqdata.toml"
	params:
		script_path=SCRIPT_PATH
	script:
		"combine_seqdata.py"

rule snippy:
    input:
        seqdata = '{sample}/seqdata.toml',
        kraken = '{sample}/kraken.toml'
        
    output:
        '{sample}/snippy.toml',
    threads:
        8
    # singularity: SINGULARITY_PATH_SNIPPY
    params:
        script_path=SCRIPT_PATH,
        reference = REFERENCE
    script:
        "snippy.py"

rule qc_snippy: 
    input:
        '{sample}/snippy.toml'
        
    output:
        '{sample}/snippy_qc.toml'
    params:
        script_path = SCRIPT_PATH,
        minaln = MIN_ALN
    script:
        "snippy_qc.py"

rule run_snippy_core:
    input:
        expand("{sample}/snippy_qc.toml", sample = SAMPLE)
    output:
        'snippy_core.toml'
    # singularity: SINGULARITY_PATH_SNIPPY
    params:
        mask_string = MASK,
        script_path = SCRIPT_PATH,
        reference = REFERENCE
    script:
        "snippy_core.py"

rule run_snpdists:
    input:
        'snippy_core.toml'
    output:
        'distances.toml' 
    # singularity:SINGULARITY_PATH_SNIPPY
    params:
        script_path = SCRIPT_PATH
    script:
        "snp_dists.py"

rule run_iqtree_core:
    input:
        core = 'snippy_core.toml', 
        ref = REFERENCE, 
        idx = IDX
    
    output:
        'iqtree.toml',
    params:
        script_path = SCRIPT_PATH	
    script:
        "run_iqtree.py"
	
rule run_tbprofiler:
    input:
        snippy = "{sample}/snippy.toml",
        kraken = "{sample}/kraken.toml",
        qc = "{sample}/snippy_qc.toml"
    # singularity:SINGULARITY_PATH_PROFILER
    params:
        script_path = SCRIPT_PATH,
        threads = PROFILER_THREADS
    output:
        "{sample}/tbprofiler.toml"
    script:
        "run_tbprofiler.py"

rule run_snpit:
    input:
        tbprofiler = "{sample}/tbprofiler.toml"
    output:
        "{sample}/snpit.toml"
    singularity:"docker://mduphl/snpit:v1.0.0"
    params:
        script_path = SCRIPT_PATH,
    script:
        "run_snpit.py"
        
rule collate_resistance:
    input:
        expand("{sample}/snpit.toml", sample = SAMPLE)
    output:
        'resistance.toml'
    params:
        script_path = SCRIPT_PATH,
        db_version = DB_VERSION,
        mode = MODE
    script:
        "collate.py"

rule collate_report:
    input:
        'resistance.toml','iqtree.toml'
    output:
        'report.toml', 'index.html'
    params:
        script_path = SCRIPT_PATH,
        template_path = TEMPLATE_PATH,
        wkdir = WORKDIR,
        amr_only = AMR_ONLY,
        idx = IDX
    script:
        "write_report.py"