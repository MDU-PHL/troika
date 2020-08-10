#!/usr/bin/env python3.7

import subprocess, pathlib, json, sys, toml

from snakemake import shell


def run_snpit(vcf):
    # print(vcf)   
    v = pathlib.Path(vcf)
    if v.exists():
        cmd = f"snpit --input {v}"
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
    # print(j)
    # d = j
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
    # print(data)
    with open(output, 'w') as f:
        
        toml.dump(data, f)
    
def main(tbprofiler,isolate):
    
    tml = open_toml(tbprofiler)
    # print(tml)
    data = {}
    data[isolate] = {}
    data[isolate]['snpit'] = {}
    data[isolate]['tb_profiler'] = {}
    if tml[isolate]['tb_profiler']['done'] == 'Yes':
        vcf = f"{isolate}/snps.raw.vcf"
        snpit = run_snpit(vcf = vcf)
        # print(snpit[0], snpit[1],snpit[1])
        data[isolate]['tb_profiler']['done'] = 'Yes'
        data[isolate]['tb_profiler']['kmer-id'] = tml[isolate]['tb_profiler']['kmer-id']
        # print(data[isolate])
        data[isolate]['snpit']['species'] = snpit[0]
        # # print(data[isolate]['snpit'])
        data[isolate]['snpit']['lineage'] = snpit[1]
        # # print(data[isolate]['snpit'])
        data[isolate]['snpit']['family'] = snpit[2]
        # print(data[isolate]['snpit'])
        data[isolate]['snpit']['done'] = 'Yes'
        # print(data[isolate]['snpit'])
        append_to_json(lineage = snpit, json_file = f'{isolate}/tbprofiler.results.json', outfile = f'{isolate}/tbprofiler.snpit_results.json')
        # print(isolate)
        
    else:
        data[isolate]['tb_profiler']['done'] = 'No'
        data[isolate]['snpit']['done'] = 'No'
        data[isolate]['tb_profiler']['kmer-id'] = tml[isolate]['tbprofiler']['kmer-id']
        data[isolate]['snpit']['data'] = {}
    write_toml(data = data, output = f"{isolate}/snpit.toml")

    
tbprofiler = snakemake.input.tbprofiler
isolate = snakemake.wildcards.sample

main(tbprofiler = tbprofiler, isolate = isolate)