import toml, pathlib, subprocess, sys, pandas

from snakemake import shell

def get_top_3(isolate):
    
    report = pathlib.Path(isolate, 'kraken2.tab')
    df= pandas.read_csv(report, sep = "\t", header =None, names = ['percentage', 'frag1', 'frag2','code','taxon','name'])
    df['percentage'] = df['percentage'].apply(lambda x:float(x.strip('%')) if isinstance(x, str) == True else float(x)) #remove % from columns
    df = df.sort_values(by = ['percentage'], ascending = False)
    df = df[df['code'].isin(['U','S'])]     
    df = df.reset_index(drop = True) 
    top_three = [
        (f"{df.loc[0,'name'].strip()}",round(df.loc[0,'percentage'],2)),
        (f"{df.loc[1,'name'].strip()}",round(df.loc[1,'percentage'],2)),
        (f"{df.loc[2,'name'].strip()}",round(df.loc[2,'percentage'],2))
    ]

    return top_three

def generate_cmd(r1, r2, isolate, db, data):
    
    
    
    cmd = f"kraken2 --paired {r1} {r2} --minimum-base-quality 13 --report {isolate}/kraken2.tab --memory-mapping --db {db} --threads 4"
    data[isolate]['kraken']['kraken_db'] = f"{db}"

    return cmd, data

def run_cmd(cmd):
    
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode

def extract_metrics(mash_string):
    
    for m in mash_string:
        if 'Estimated coverage' in m:
            d = m.split(':')[-1].strip()
    return d
    

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(r1, r2, isolate, kraken_db,run_kraken):
    
    # set up data dict
    tml = f"{pathlib.Path(isolate, 'kraken.toml')}"
    data = {}
    data[isolate] = {}
    data[isolate]['kraken'] = {}
    # run kraken
    if run_kraken:
        cmd, data = generate_cmd(r1 = r1, r2 = r2, isolate = isolate, db = kraken_db, data = data)
        kraken_returncode = run_cmd(cmd)
        if kraken_returncode == 0:
            # add to data dict
            data[isolate]['kraken']['done'] = 'Yes'
            data[isolate]['kraken']['file'] = f"{isolate}/kraken2.tab"
            data[isolate]['kraken']['top_3'] = get_top_3(isolate)
    else:
        data[isolate]['kraken']['done'] = 'No'
    write_toml(data = data, output= f"{isolate}/kraken.toml")

r1 = snakemake.input.r1
r2 = snakemake.input.r2
isolate = snakemake.wildcards.sample
run_kraken = snakemake.params.run_kraken
kraken_db = snakemake.params.kraken_db

main(r1 = r1, r2 = r2, isolate = isolate,run_kraken = run_kraken, kraken_db = kraken_db)