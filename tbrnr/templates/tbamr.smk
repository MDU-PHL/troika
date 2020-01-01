configfile:'config.yaml'

SAMPLES = config['samples'].split()
SAMPLE_STRING = ','.join(SAMPLES)
SINGULARITY_PATH = config['singularity_path'] # path to container for running tb-profiler
THREADS = config['threads']
SCRIPT_PATH = config['script_path']
FINAL_OUTPUT = config['final_output'] # a two item list - jobid.csv, jobid.json
DB_VERSION = config['db_version']

rule all:
    input:
        FINAL_OUTPUT

rule run_tbprofiler:
    input:
        r1 = "READS/{sample}/R1.fq.gz",
        r2 = "READS/{sample}/R2.fq.gz"
    singularity:SINGULARITY_PATH
    threads:THREADS
    output:
        "results/{sample}.results.json", "vcf/{sample}.targets.vcf.gz"
    shell:
        """
        tb-profiler profile -1 {input.r1} -2 {input.r2} -p {wildcards.sample} --caller BCFtools -t {threads}
        """

rule collate_resistance:
    input:
        expand("results/{sample}.results.json", sample = SAMPLES)
    output:
        'tbrnr.csv', 'tbrnr.jsob'
    params:
        script_path = SCRIPT_PATH
        db_version = DB_VERSION
        sample = SAMPLE_STRING
    shell:
        """
        python3 {params.script_path}/collate.py {params.sample} {params.db_version}
        """

rule mv_outputs:
    input:
        'tbrnr.csv', 'tbrnr.json'
    output:
        FINAL_OUTPUT
    shell:
        """
        mv {input[0]} {output[0]}
        mv {input[1]} {output[1]}
        """