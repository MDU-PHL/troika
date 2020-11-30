import toml, pathlib, subprocess, sys
from snakemake import shell


def generate_snippy_cmd(r1, r2, isolate, reference, threads):
    
    cmd = f"snippy --outdir {isolate} --ref {reference} --R1 {r1} --R2 {r2} --force --cpus {threads}"
    print(cmd)
    return cmd

def run_cmd(cmd):
    
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode

def get_species(kraken):
    mtbc = ['pinnipedii', 'tuberculosis', 'orygis', 'bovis', 'bovis BCG', 'canetti', 'microti', 'africanum']
    tml = open_toml(kraken)
    print(tml)
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

def get_reads(seqdata, isolate):

    s = open_toml(seqdata)
    r1 = s[isolate]['seqdata']['R1']
    r2 = s[isolate]['seqdata']['R2']
    return r1,r2

def get_quality(seqdata, kraken, isolate):

    s = open_toml(seqdata)
    print(s)
    k = get_species(kraken)
    print(k)
    if s[isolate]['seqdata']['data']['Quality'] == 'PASS' and k == 'MTBC':
        return 'Yes'
    else:
        return 'No'


def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(seqdata,kraken, isolate, output, reference, threads):
    
    print(seqdata)
    # print(threads)
    r1,r2 = get_reads(seqdata = seqdata, isolate = isolate)
    run_snippy = get_quality(seqdata = seqdata, kraken = kraken, isolate = isolate)
    print(run_snippy)
    data = {}
    data[isolate] = {}
    data[isolate]['snippy'] = {}
    data[isolate]['snippy']['reference'] = reference
    data[isolate]['snippy']['R1'] = r1
    data[isolate]['snippy']['R2'] = r2
    data[isolate]['snippy']['run_snippy'] = run_snippy

    if run_snippy == 'Yes':
        cmd = generate_snippy_cmd(r1 = r1, r2=r2, isolate = isolate, reference = reference, threads = threads)
        # print(cmd)
        p = run_cmd(cmd)
        # print(p)
        if p == 0:
            data[isolate]['snippy']['bam'] = f"{isolate}/snps.bam"
            data[isolate]['snippy']['alignment'] = f"{isolate}/snps.aligned.fa"
            data[isolate]['snippy']['vcf'] = f"{isolate}/snps.vcf"
            data[isolate]['snippy']['cmd'] = cmd
    write_toml(data = data, output = f"{isolate}/snippy.toml") 

seqdata = snakemake.input.seqdata
kraken = snakemake.input.kraken
isolate = snakemake.wildcards.sample
output = snakemake.output
reference = snakemake.params.reference
threads = snakemake.threads

main(seqdata = seqdata, kraken =kraken, isolate = isolate, output = output,reference = reference, threads =threads)
