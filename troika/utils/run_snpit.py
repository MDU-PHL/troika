import subprocess, pathlib, json, sys

# print(vcfs)
def unzip_vcf(vcf):
    vcf_id = vcf.name.split('.')[0]
    if pathlib.Path(vcf).suffix == ".gz":
        gunzip= f"gunzip -d -c {vcf} > vcf/{vcf_id}.vcf"
        subprocess.run(gunzip, shell = True, capture_output = True)
        return f"vcf/{vcf_id}.vcf"
    else:
        return vcf

def run_snpit(vcf):
    # print(vcf)   
    v = pathlib.Path(vcf)
    if v.exists():
        unzipped_vcf = unzip_vcf(vcf = v)
        cmd = f"/home/khhor/.local/bin/snpit --input {unzipped_vcf}" # convert to non kristy specific
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



if __name__ == "__main__":
    # print(len(sys.argv))
    if len(sys.argv) == 4:
        vcf = sys.argv[1]
        # print(vcf)
        json_file = sys.argv[2]
        outfile = sys.argv[3]
        lineage = run_snpit(vcf)
        append_to_json(lineage, json_file, outfile)

