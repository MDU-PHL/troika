import pandas, pathlib, subprocess

configfile:'config.yaml'

SAMPLE = config['samples'].split()
SINGULARITY_PATH_PROFILER = config['singularity_path_profiler'] # path to container for running tb-profiler
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
        'kraken.toml', 'iqtree.toml', 'report.toml', 'seqdata.toml', 'distances.toml','snippy_core.toml','resistance.toml', expand("{sample}/snpit.toml", sample = SAMPLE)
        

# if RUN_SPECIES == True:
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
    shell:
        """
        python3 {params.script_path}/kraken.py {input.r1} {input.r2} {wildcards.sample} {params.run_kraken} {params.kraken_db} 
        """
rule combine_kraken:
    input:
        expand("{sample}/kraken.toml", sample = SAMPLE)
    output:
        "kraken.toml"
    params:
        script_path = SCRIPT_PATH
    shell:
        """
        python3 {params.script_path}/combine_kraken.py {input}
        """

rule estimate_coverage:
	input:
		r1="{sample}/R1.fq.gz",
		r2="{sample}/R2.fq.gz"
	output:
		"{sample}/mash.toml"
	params:
		script_path = SCRIPT_PATH
	
	shell:
		"""
		python3 {params.script_path}/mash.py {input.r1} {input.r2} {wildcards.sample} {output}
		"""

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
	
	shell:
		"""
		python3 {params.script_path}/seqdata.py {input.r1} {input.r2} {wildcards.sample} {input.mash} {params.mincov}
		"""

rule combine_seqdata:
	input:
		expand("{sample}/seqdata.toml", sample = SAMPLE)
	output:
		"seqdata.toml"
	params:
		script_path=SCRIPT_PATH
	shell:
		"""
		python3 {params.script_path}/combine_seqdata.py {input} 
		"""
rule snippy:
    input:
        seqdata = '{sample}/seqdata.toml',
        kraken = '{sample}/kraken.toml'
        
    output:
        '{sample}/snippy.toml',
    threads:
        8
    
    params:
        script_path=SCRIPT_PATH,
        reference = REFERENCE
    shell:
        """
        python3 {params.script_path}/snippy.py {input.seqdata} {input.kraken} {wildcards.sample} {output} {params.reference} {threads}
        """
    

rule qc_snippy: 
    input:
        '{sample}/snippy.toml'
        
    output:
        '{sample}/snippy_qc.toml'
    params:
        script_path = SCRIPT_PATH,
        minaln = MIN_ALN
    shell:
        """
        python3 {params.script_path}/snippy_qc.py {input} {wildcards.sample} {output} {params.minaln}
        """

rule run_snippy_core:
    input:
        expand("{sample}/snippy_qc.toml", sample = SAMPLE)
    output:
        'snippy_core.toml'
    
    params:
        mask_string = MASK,
        script_path = SCRIPT_PATH,
        reference = REFERENCE
    shell:
        """
        python3 {params.script_path}/snippy_core.py {params.mask_string} {params.reference} {input}
        """

rule run_snpdists:
    input:
        'snippy_core.toml'
    output:
        'distances.toml' 
    params:
        script_path = SCRIPT_PATH
    shell:
        """
        python3 {params.script_path}/snp_dists.py {input}
        """

rule run_iqtree_core:
    input:
        core = 'snippy_core.toml', 
        ref = REFERENCE, 
        idx = IDX
    
    output:
        'iqtree.toml',
    params:
        script_path = SCRIPT_PATH	
    shell:
        """	
        python3 {params.script_path}/run_iqtree.py {input.core} {input.ref} {input.idx} {params.script_path}
        """
	
rule run_tbprofiler:
    input:
        snippy = "{sample}/snippy.toml",
        kraken = "{sample}/kraken.toml",
        qc = "{sample}/snippy_qc.toml"
    singularity:SINGULARITY_PATH_PROFILER
    params:
        script_path = SCRIPT_PATH,
        threads = PROFILER_THREADS
    output:
        "{sample}/tbprofiler.toml"
    shell:
        """
        python3 {params.script_path}/run_tbprofiler.py {input.snippy} {input.qc} {input.kraken} {wildcards.sample} {params.threads}
        """

rule run_snpit:
    input:
        # snippy = "{sample}/snippy.toml",
        tbprofiler = "{sample}/tbprofiler.toml"
    output:
        "{sample}/snpit.toml"
    params:
        script_path = SCRIPT_PATH,
    shell:
        """
        python3 {params.script_path}/run_snpit.py {input.tbprofiler} {wildcards.sample} 
        """
        
rule collate_resistance:
    input:
        expand("{sample}/snpit.toml", sample = SAMPLE)
    output:
        'resistance.toml'
    params:
        script_path = SCRIPT_PATH,
        db_version = DB_VERSION,
        mode = MODE
    shell:
        """
        python3 {params.script_path}/collate.py {params.db_version} {params.mode} {input}
        """

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
    shell:
        """
        python3 {params.script_path}/write_report.py {params.template_path} {params.wkdir} {params.amr_only} {params.idx}
        """