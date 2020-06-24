import toml, pathlib, subprocess, sys, pandas, datetime, json

from snakemake import shell

def get_species(kraken):
    mtbc = ['pinnipedii', 'tuberculosis', 'orygis', 'bovis', 'bovis BCG', 'canetti', 'microti', 'africanum']
    tml = open_toml(kraken)
    # print(tml)
    isolate = list(tml.keys())[0]
    if tml[isolate]['kraken']['done'] == 'No':
        return 'kmer-id not performed'
    elif tml[isolate]['kraken']['done'] == 'Yes':
        s = tml[isolate]['kraken']['top_3'][0][0]
        # print(s)
        if 'Mycobacterium' in s:
            print('Mycobacterium genus found')
            for m in mtbc:
                print(m)
                print(m in s)
                if m in s:
                    print("found MTBC species")
                    return 'MTBC'
    # return 'non-MTBC'
            

def generate_tbprofiler_cmd(snippy, isolate, threads):
    tml = open_toml(snippy)
    bam = f"{tml[isolate]['snippy']['bam']}"
    cmd = f"tb-profiler profile -a {bam} --caller bcftools -t 6  -d {isolate}"

    return cmd

def generate_mv_cmd(isolate):

    cmd = f"mv {isolate}/results/tbprofiler.results.json {isolate}/tbprofiler.results.json"
    return cmd

def get_data(isolate):
    
    if pathlib.Path(isolate,"tbprofiler.results.json").exists():
        d = json.load(open(f"{isolate}/tbprofiler.results.json"))
        return d
    else:
        return {}

def generate_rm_cmd(isolate):

    cmd = f"rm -rf {isolate}/results {isolate}/bam {isolate}/vcf"
    return cmd

def run_cmd(cmd):

    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(snippy,qc, isolate, threads, kraken):
    
    # print(isolate)
    
    species = get_species(kraken)
    tml = open_toml(snippy)
    # print(tml)
    data = {}
    data[isolate] = {}
    data[isolate]['tbprofiler'] = {}
    
    if tml[isolate]['snippy']['run_snippy'] == 'Yes' and species in ['MTBC', 'kmer-id not performed']:
        cmd = generate_tbprofiler_cmd(snippy = snippy, isolate = isolate, threads = threads)
        print(cmd)
        p = run_cmd(cmd)
        if p == 0:
            run_cmd(generate_mv_cmd(isolate = isolate))
            # run_cmd(generate_rm_cmd)
            data[isolate]['tbprofiler']['done'] = 'Yes'
            data[isolate]['tbprofiler']['kmer-id'] = species
            data[isolate]['tbprofiler']['data'] = get_data(isolate=isolate)
    else:
        data[isolate]['tbprofiler']['done'] = 'No'
        data[isolate]['tbprofiler']['kmer-id'] = 'kmer id not consistent with MTBC species'
        data[isolate]['tbprofiler']['data'] = get_data(isolate=isolate) 
    write_toml(data = data, output = f"{isolate}/tbprofiler.toml")   

snippy = snakemake.input.snippy
qc = snakemake.input.qc
kraken = snakemake.input.kraken
isolate = snakemake.wildcards.sample
threads = snakemake.params.threads

main(snippy = snippy, qc = qc, kraken = kraken, isolate = isolate, threads = threads)
