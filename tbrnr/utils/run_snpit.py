import subprocess, pathlib, json

# print(vcfs)
def unzip_vcf(vcf):
    vcf_id = vcf.name.split('.')[0]
    if pathlib.Path(vcf).suffix == ".gz":
        gunzip= f"gunzip -d -c {vcf} > vcf/{vcf_id}.vcf"
        subprocess.run(gunzip, shell = True, capture_output = True)
        return f"vcf/{vcf_id}.vcf"

def run_snpit(vcf):   
    v = pathlib.Path(vcf)

    if v.exists():
        unzipped_vcf = unzip_vcf(vcf = v)
        cmd = f"snpit-run.py --input {unzipped_vcf}"
        p = subprocess.run(cmd, shell = True, capture_output = True, encoding = "utf-8")
        # print(p.stdout)
        lineage_strs = p.stdout.split("    ")
        # print(lineage_strs)
        if not "tuberculosis" in lineage_strs[0]:
            sp = "Non tuberculosis MTBC"
            l = lineage_strs[1]
            perc = lineage_strs[-1].strip()
        elif "tuberculosis" in lineage_strs[0]:
            sp = "M. tuberculosis"
            l = [i for i in lineage_strs if "Lineage" in i][0]
            perc = lineage_strs[-1].strip()
        return [sp, l]

def append_to_json(lineage, json_file):
    isolate = json_file.name.split('.')[0]
    j = json.load(open(json_file))
    j['snpit'] = {
        'species':lineage[0],
        'lineage':lineage[1]
    }
    json.dump(j, f"results/{isolate}.snpit_results.json")



if __name__ == "__main__":

    if len(sys.argv) == 3:
        vcf = sys.argv[1]
        json_file = sys.argv[2]
        lineage = run_snpit(vcf)
        append_to_json(lineage, json_file)

