import pandas, pathlib, subprocess

configfile:'config.yaml'


def get_report_input(amr_only, final_output):
    
    fn = final_output.split(',')
    fn = [f.strip("''\" ") for f in fn]
    # print(fn)
    report_input = [f"'report/{f.strip()}'" for f in fn if 'core.txt' not in f]

    if not amr_only:
        report_input.append(f"'report/core_genome.tab'")
    
    return ','.join(report_input)

def get_final_output_string(string, amr_only):
    string = ','.join([f"'{s}'" for s in string.split(',')])
    file_list = f"'seqdata.tab',{string}"
    if amr_only:
        return file_list
    else:
        file_list = f"{file_list}, 'core.treefile', 'distances.tab','core.tab', 'core.txt'"
        return file_list
def get_troika_output(string):
    s = string.split(',')
    troika = ','.join([f"'{t}'" for t in s if 'mtb.tab' not in t])

    return troika

def get_troika_report(string):
    sn = string.split(',')
    sn = [s.strip("''\" ") for s in sn]
    report_string = ','.join([f"'report/{s}'" for s in sn])
    return report_string


SAMPLES = config['samples'].split()
SAMPLE_STRING = ','.join(SAMPLES)
print(SAMPLE_STRING)
SAMPLE_FOR_CORE = ' '.join(SAMPLES)
SINGULARITY_PATH_PROFILER = config['singularity_path_profiler'] # path to container for running tb-profiler
PROFILER_THREADS = config['profiler_threads']
SCRIPT_PATH = config['script_path']
 # a two (or three item) item list - jobid.csv, jobid.json, mtb.tab

DB_VERSION = config['db_version']
RUN_SPECIES = config['run_species']
KRAKEN_DB = config['kraken_db']
# PREFILL_PATH = config['prefill_path']
KRAKEN_THREADS = config['kraken_threads']
SNIPPY_THREADS = config['snippy_threads']
REFERENCE = config['reference']
MASK = config['mask']
AMR_ONLY = config['amr_only']
WORKDIR = f"'{pathlib.Path.cwd().absolute()}'"
TEMPLATE_PATH = config['template_path']
TROIKA_OUTPUT = get_troika_output(config['final_output'])
TROIKA_REPORT = get_troika_report(string= TROIKA_OUTPUT)
FINAL_OUTPUT = get_final_output_string(config['final_output'], amr_only = AMR_ONLY)
REPORT_INPUT = get_report_input(amr_only = AMR_ONLY, final_output= FINAL_OUTPUT)
print(REPORT_INPUT)
print(TROIKA_OUTPUT)
print(TROIKA_REPORT)

rule all:
    input:
        # REPORT_INPUT,
        'report/seqdata.tab','report/troika.tab','report/troika.json','report/core.treefile','report/distances.tab','report/core.tab','report/core_genome.tab',
        expand("{sample}/tbprofiler.snpit_results.json", sample = SAMPLES), 
        'core.aln',
        'report/index.html'


if RUN_SPECIES == True:
    rule species:
        input:
            r1 = "{sample}/R1.fq.gz",
            r2 = "{sample}/R2.fq.gz",
        output:
            '{sample}/kraken_report.txt'
        params:
            db = KRAKEN_DB,
            threads = KRAKEN_THREADS
        shell:
            """
            kraken2 --paired {input.r1} {input.r2} --db {params.db} --threads {params.threads} --minimum-base-quality 13 --report {output} --classified-out /dev/null --unclassified-out /dev/null
            """

    rule combine_kraken:
        input:
            expand("{sample}/kraken_report.txt", sample = SAMPLES)
        output:
            "mtb.tab"
        params:
            script_path = SCRIPT_PATH
        shell:
            """
            python3 {params.script_path}/combine_kraken.py {input} {output}
            """
            
rule seqdata:
    input:
        '{sample}/R1.fq.gz',
        '{sample}/R2.fq.gz'
    output:
        "{sample}/seqdata.tab"
    shell:
        """
        seqtk fqchk {input[0]} {input[1]} > {output}
        """

rule estimate_coverage:
    input:
        "{sample}/R1.fq.gz",
        "{sample}/R2.fq.gz"
    output:
        "{sample}/mash.txt"
    shell:
        """
        mash sketch -r {input[0]} {input[1]} -m 3 -k 31 -o mash  &> {output}
        """

rule generate_yield:
    input:
        "{sample}/mash.txt",
        "{sample}/seqdata.tab"
    output:
        "{sample}/yield.tab"
    params:
        script_path = SCRIPT_PATH
    shell:
        """
        python3 {params.script_path}/generate_yield.py {input[1]} {input[0]} {output}
        """

rule combine_seqdata:
    input:
        expand("{sample}/yield.tab", sample = SAMPLES)
    output:
        "seqdata.tab"
    run:
        import pathlib, pandas, numpy
        sdfiles = f"{input}".split()
        seqdata = pandas.DataFrame()
        for sd in sdfiles:
            p = pathlib.Path(sd)
            df = pandas.read_csv(sd, sep = "\t")
            df['Isolate'] = f"{p.parts[0]}"
            if seqdata.empty:
                seqdata = df
            else:
                seqdata = seqdata.append(df)
        seqdata['Quality'] = numpy.where(seqdata['Estimated depth'] >= 40, 'PASS','FAIL')
        seqdata = seqdata[['Isolate','Reads','Yield','GC content','Min len','Avg len','Max len','Avg Qual','Estimated depth', 'Quality']]
        seqdata.to_csv(f"{output}", sep = '\t', index = False)
    
rule run_snippy:
    input:
        r1 = "{sample}/R1.fq.gz",
        r2 = "{sample}/R2.fq.gz",
    output:
        bam = '{sample}/snps.bam',
        raw_vcf = '{sample}/snps.raw.vcf',
        aln = '{sample}/snps.aligned.fa'
    params:
        threads = SNIPPY_THREADS,
        reference = REFERENCE
    
    shell:
        """
        snippy --outdir {wildcards.sample} --ref {params.reference} --R1 {input.r1} --R2 {input.r2} --force --cpus {params.threads}
        """

# if not AMR_ONLY:

rule run_snippy_core:
    input:
        expand('{sample}/snps.aligned.fa', sample = SAMPLES)
    output:
        'core.vcf',
        'core.txt',
        'core.aln', 
        'core.full.aln',
        'core.tab'
    params:
        reference = REFERENCE,
        sample_list = SAMPLE_FOR_CORE,
        mask = MASK
    shell:
        """
        snippy-core --ref {params.reference} --mask {params.mask} {params.sample_list}
        """

rule run_snpdists:
    input:
        'core.aln'
    output:
        'distances.tab' 
    shell:
        """
        snp-dists {input} > {output}
        """
        
rule index_reference:
    input:
        REFERENCE
    output:
        "ref.fa",
        "ref.fa.fai"
    run:
        from Bio import SeqIO
        import pathlib, subprocess
        ref = f"{output[0]}"
        idx = f"{output[1]}"
        if '.fa' not in REFERENCE:
            print(f"converting {REFERENCE}")
            SeqIO.convert(f"{input[0]}", 'genbank', ref    , 'fasta')
            print(f"converted {REFERENCE}")
        else:
            subprocess.run(f"ln -sf {REFERENCE} {ref}", shell = True)
        subprocess.run(f"samtools faidx {ref}", shell =True)


rule calculate_iqtree_command_core:
    input:
        aln = 'core.aln',
        ref = "ref.fa"
    output:
        'run_iqtree_core.sh'
    params:
        script_path = SCRIPT_PATH,
    shell:
        "bash {params.script_path}/iqtree_generator.sh {input.ref} {input.aln} core 20 > {output}"



rule run_iqtree_core:
    input:
        'run_iqtree_core.sh'
    
    output:
        'core.iqtree',
        'core.treefile',
    shell:
        """    
        bash run_iqtree_core.sh
        rm *.ckp.gz *.contree *.bionj
        """ 


rule collate_report_files:
    input:
        'seqdata.tab', 'core.txt', 'core.treefile', 'core.tab', 'distances.tab', 'core.tab'
    output:
        'report/seqdata.tab', 'report/core_genome.tab', 'report/core.treefile','report/distances.tab','report/core.tab'
    run:        
        import pandas, pathlib, subprocess, numpy
        
        # for core.txt
        df = pandas.read_csv(pathlib.Path(f"core.txt"), sep = '\t')
        df['% USED'] = 100 * (df['LENGTH'] - df['UNALIGNED'])/ df['LENGTH']
        df['% USED'] = df['% USED'].round(2)
        df = df.rename(columns={'ID':'Isolate'})
        df.to_csv(f"report/core_genome.tab", sep='\t', index = False)

        cmd = f"""
cp seqdata.tab report/seqdata.tab
cp core.treefile report/core.treefile
cp distances.tab report/distances.tab
cp core.tab report/core.tab
"""
        subprocess.run(cmd, shell = True)

rule run_tbprofiler:
    input:
        bam = "{sample}/snps.bam"
    singularity:SINGULARITY_PATH_PROFILER
    threads:PROFILER_THREADS
    output:
        "{sample}/tbprofiler.results.json"
    shell:
        """
        tb-profiler profile -a {input.bam} --caller bcftools -t {threads} -d {wildcards.sample}
        mv {wildcards.sample}/results/tbprofiler.results.json {output}
        rm -rf {wildcards.sample}/results {wildcards.sample}/bam {wildcards.sample}/vcf
        """

rule run_snpit:
    input:
        vcf = "{sample}/snps.raw.vcf",
        json = "{sample}/tbprofiler.results.json"
    output:
        "{sample}/tbprofiler.snpit_results.json"
    params:
        script_path = SCRIPT_PATH,
    shell:
        """
        python3 {params.script_path}/run_snpit.py {input.vcf} {input.json} {output}
        """
        
rule collate_resistance:
    input:
        expand("{sample}/tbprofiler.snpit_results.json", sample = SAMPLES)
    output:
        'troika.tab', 'troika.json'
    params:
        script_path = SCRIPT_PATH,
        db_version = DB_VERSION,
        sample = SAMPLE_STRING
    shell:
        """
        python3 {params.script_path}/collate.py {params.sample} {params.db_version}
        """

rule mv_outputs:
    input:
        'troika.tab', 'troika.json'
    output:
        'report/troika.tab', 'report/troika.json'
    shell:
        """
        mv {input[0]} {output[0]}
        mv {input[1]} {output[1]}
        """

rule collate_report:
    input:
        'report/seqdata.tab','report/troika.tab','report/troika.json','report/core.treefile','report/distances.tab','report/core.tab','report/core_genome.tab',
    output:
        'report/index.html'
    params:
        script_path = SCRIPT_PATH,
        template_path = TEMPLATE_PATH,
        wkdir = WORKDIR,
        amr_only = AMR_ONLY
    shell:
        """
        python3 {params.script_path}/write_report.py {params.template_path} {params.wkdir} {params.amr_only}
        """