import subprocess, pathlib, json, sys, toml

from snakemake import shell
def run_snpit(vcf):
    # print(vcf)   
    v = pathlib.Path(vcf)
    if v.exists():
        cmd = f"/home/khhor/.local/bin/snpit --input {v}" # convert to non kristy specific
        p = subprocess.run(cmd, shell = True, capture_output = True, encoding = "utf-8")
        snpit_data = p.stdout.split("\n")
        lineage_strs = snpit_data[1].split('\t')
        sp = lineage_strs[1]
        l = lineage_strs[2]
        f = lineage_strs[4] if lineage_strs[4] != 'N/A' else ''
        
        return [sp, l,f]

def append_to_json(lineage, json_file, outfile):
    isolate = pathlib.Path(json_file).name.split('.')[0]
    j = json.load(open(json_file))
    j['snpit'] = {
        'species':lineage[0],
        'lineage':lineage[1],
        'family':lineage[2]
    }

    with open(outfile, 'w') as f:
        json.dump(j, f)
    return j

def run_cmd(cmd):

    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(tbprofiler,isolate):
    
    tml = open_toml(tbprofiler)
    print(tml)
    data = {}
    data[isolate] = {}
    data[isolate]['snpit'] = {}
    data[isolate]['tbprofiler'] = {}
    if tml[isolate]['tbprofiler']['done'] == 'Yes':
        vcf = f"{isolate}/snps.raw.vcf"
        snpit = run_snpit(vcf = vcf)
        data[isolate]['tbprofiler']['done'] = 'Yes'
        data[isolate]['tbprofiler']['kmer-id'] = tml[isolate]['tbprofiler']['kmer-id']
        data[isolate]['snpit']['species'] = snpit[0]
        data[isolate]['snpit']['lineage'] = snpit[1]
        data[isolate]['snpit']['family'] = snpit[2]
        data[isolate]['snpit']['done'] = 'Yes'
        data[isolate]['snpit']['data'] = append_to_json(lineage = snpit, json_file = f'{isolate}/tbprofiler.results.json', outfile = f'{isolate}/tbprofiler.snpit_results.json')
    else:
        data[isolate]['tbprofiler']['done'] = 'No'
        data[isolate]['snpit']['done'] = 'No'
        data[isolate]['tbprofiler']['kmer-id'] = tml[isolate]['tbprofiler']['kmer-id']
        data[isolate]['snpit']['data'] = {}
    write_toml(data = data, output = f"{isolate}/snpit.toml")

    
tbprofiler = snakemake.input.tbprofiler
isolate = snakemake.wildcards.sample

main(tbprofiler = tbprofiler, isolate = isolate)