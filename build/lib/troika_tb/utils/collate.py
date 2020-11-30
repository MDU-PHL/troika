import json,pathlib,pandas, sys, datetime, toml, numpy

from snakemake import shell

def return_species(data):
    '''
    data is the dictionary from json.load
    '''
    
    if data['snpit']['species'] != 'N/A':
        return data['snpit']['species']
    elif 'BOV' not in data['lineage'][0]['lin']:
        return 'M. tuberculosis'
    elif 'BOV' in data['lineage'][0]['lin']:
        return 'non-tuberculosis MTBC'
    else:
        return "species undefined"
    

def return_lineage(data):
    '''
    data is the dictionary from json.load
    '''

    if data['lineage'] != []:
        main_lineage = data['main_lin'].replace('lineage', 'Lineage ')
        semicolon_check = main_lineage.split(';')
        if len(semicolon_check) > 1 :
            lins = [s.split('.')[0] for s in semicolon_check]
            if len(set(lins)) > 1:
                lineage_string = 'phylogenetic lineage undefined' 
            else:
                f = data['lineage'][0]['family'].split()[0]
                lineage_string = f"{f} ({lins[0]})"
        else:
            lineage_string = f"{data['lineage'][0]['family']} ({main_lineage})"
    elif data['lineage'] == [] and data['snpit']['lineage'] != 'N/A':
        lineage_string = f"{data['snpit']['family']} ({data['snpit']['lineage']})"
    else:
        lineage_string = 'phylogenetic lineage undefined' 

    return lineage_string

def construct_string(mutation):

    gene = mutation['gene_name'] if 'gene_name' in mutation else mutation['gene']
    change = mutation['change']
    gene = f"{gene[0].upper()}{gene[1:]}"
    return f"{gene} ({change})"    

def check_variants(mutation, drug):
    '''
    mutation is the dictionary from dr_variants
    based on validation of mms155
    allelic frequency minimum set to .1 so no need to filter
    for Rif, Inh and Emb - use mutations that are NOT indeterminate
    for all other drugs - use ALL mutations in the dr_variants - ie include the indeterminate
    '''
    
    exclude_ind = ['rifampicin', 'isoniazid', 'ethambutol']
    if drug in exclude_ind and mutation['confidence'] != 'indeterminate':
        return construct_string(mutation = mutation)
    elif drug not in exclude_ind:
        return construct_string(mutation = mutation)
    else:
        return ''

def get_dr_variants(data):
    '''
    data is the dictionary from json.load
    if drug in 'rifampicin', 'isoniazid', 'ethambutol' the EXCLUDE 'indeterminate' mutations, eles include them 
    '''
    # Check with Norelle - whether or not she wants to include drugs which did not fall in the validated dataset
    # a dictionary defaulted to 'No mutation detected'
    d = {
        'Amikacin':'No mutation detected', #Block 2
        # 'Aminoglycosides':'No mutation detected',
        'Bedaquiline':'No mutation detected',
        'Capreomycin':'No mutation detected',
        'Ciprofloxacin':'No mutation detected', 
        'Clofazimine':'No mutation detected',
        'Cycloserine':'No mutation detected',
        'Delamanid':'No mutation detected', #Block 3
        'Ethambutol':'No mutation detected', 
        'Ethionamide':'No mutation detected', 
        'Fluoroquinolones':'No mutation detected', 
        'Isoniazid':'No mutation detected',
        'Kanamycin':'No mutation detected', 
        'Levofloxacin':'No mutation detected', 
        'Linezolid':'No mutation detected', 
        'Moxifloxacin':'No mutation detected', #Block 1
        'Ofloxacin':'No mutation detected', 
        'Para-aminosalicylic_acid':'No mutation detected',
        'Pyrazinamide':'No mutation detected',
        'Rifampicin':'No mutation detected',
        'Streptomycin': 'No mutation detected'}
    
    
    if data['dr_variants'] != []:
        for mut in data['dr_variants']:
            drug = mut['drug']
            mutation = check_variants(mutation = mut, drug = drug)
            if mutation != '' and drug not in ['aminoglycosides']: # TODO streptomycin needs to be added back in
                if d[f'{drug.capitalize()}'] != 'No mutation detected':
                    d[f'{drug.capitalize()}'] = f"{ d[f'{drug.capitalize()}']},{mutation}"
                else:
                    d[f'{drug.capitalize()}'] = mutation
    
    return d

def calculate_resistance(drugs_dict):
    '''
    based on WHO criteria for TB resistance categories
    No drug resistance predicted = no mutations in first line drugs
    Mono-resistance predicted =resistance in one of rif, inh, emb or pza
    Poly-resistance predicted  = two or more resisance where both rif and inh are not included
    Multi-drug resistance predicted = rif and inh resistance plus or minus emb and pza
    Extensive drug-resistance predicted rif and inh AND a flouroquinolone AND one of amikacin/capreomycin/kanamycin/sm #TODO add Streptomycin to logic
    '''
    # list of drigs in flq or ack for xdr
    flqlist = ['Fluoroquinolones','Levofloxacin','Moxifloxacin','Ofloxacin','Ciprofloxacin']
    ackslist = ['Amikacin', 'Capreomycin', 'Kanamycin', 'Streptomycin']
    # initialise score 
    score = 0
    # initialise flq and ack values
    flq = False
    ack = False
    # generate values for resistance call
    for d in drugs_dict: 
        if drugs_dict[d] != 'No mutation detected':
            if d in flqlist:
                flq = True
            elif d in ackslist:
                ack = True
            elif d in ['Rifampicin', 'Isoniazid']:
                score = score + 3
            elif d in ['Pyrazinamide', 'Ethambutol']:
                score = score + 1
    
    # initialise the resistance to default 
    resistance = 'No drug resistance predicted'
    if score == 1 or score == 3: # one first line drug
        resistance = 'Mono-resistance predicted'
    elif score == 2 or score in range(4,6): # more than one first line drug where INH OR RIF can be present
        resistance = 'Poly-resistance predicted'
    elif score in range(6,9): # RIF and INH +/- PZA or EMB
        if flq and ack: # there are mutations in a FLQ and amikacin,capreomycin,kanamycin
            resistance = 'Extensive drug-resistance predicted'
        else:
            resistance = 'Multi-drug resistance predicted'
    
    return resistance

def make_df(result_dict):

    cols_list = ['Isolate', 'Organism identification by WGS', 'Phylogenetic lineage', 'Predicted drug resistance','Rifampicin', 'Isoniazid', 'Pyrazinamide', 'Ethambutol','Moxifloxacin','Amikacin', 'Cycloserine', 'Ethionamide', 'Para-aminosalicylic acid','Clofazimine', 'Delamanid', 'Bedaquiline', 'Linezolid','Database version used for analysis']
    
    df = pandas.DataFrame.from_dict(data = result_dict, orient = 'index')
    df = df.reset_index()
    df = df.rename(columns = {'index': 'Isolate', 'Para-aminosalicylic_acid':'Para-aminosalicylic acid'})
    df = df.reindex(cols_list, axis = 'columns')
    df = df.fillna('No mutation detected')
    # print(df)
    return df
# ERR2120246
def save_results(result_dict):

    write_toml(data = result_dict, output = 'resistance.toml')
    df = make_df(result_dict)
    df.to_csv(f'troika.tab', index = False, sep = '\t')


def mdu_troika(result_dict):
    date = datetime.datetime.today().strftime("%d_%m_%y")
    edit_list = ['Phylogenetic lineage', 'Predicted drug resistance','Rifampicin', 'Isoniazid', 'Pyrazinamide', 'Ethambutol','Moxifloxacin','Amikacin', 'Cycloserine', 'Ethionamide', 'Para-aminosalicylic acid','Clofazimine', 'Delamanid', 'Bedaquiline', 'Linezolid']
    df = make_df(result_dict)
    for e in edit_list:
        df[e] = numpy.where(df['Organism identification by WGS'] != 'M. tuberculosis', '', df[e])
    df = df.rename(columns = {'Isolate': 'MDU sample ID',
                               'Organism identification by WGS': 'Identification (WGS)',
                               'Predicted drug resistance':'Predicted drug resist. summary:'})
    df = df.replace('No mutation detected','No mutn det')
    df.to_csv(f'MMS155_{date}.csv', index = False)

def collate_results(isolates, db, mode):
    
    # iszolates = isolates.split(',')
    result_dict = {}
    for isolate in isolates:
        p = pathlib.Path(isolate)

        json_path = p / f"tbprofiler.snpit_results.json"
        print(json_path)
        data = json.load(open(json_path))
        print(data)
        species = return_species(data)
        print(species)
        lineage = return_lineage(data)
        print(lineage)
        drugs = get_dr_variants(data)
        print(drugs)
        resistance = calculate_resistance(drugs)
        print(resistance)
        
        d = {
            'Organism identification by WGS':species,
            'Phylogenetic lineage':lineage,
            'Database version used for analysis':db,
            'Predicted drug resistance':resistance
            }
        for drug in drugs:
            d[drug] = drugs[drug]
        result_dict[isolate] = d
    # output_name = output[0].split('.')[0]
    save_results(result_dict)
    if mode == 'mdu':
        mdu_troika(result_dict)

    # if mdu:
    #     results_mdu = mdu_troika(result_dict)



# results to be recorded by tbrnr
# Species -> from snpit M. tuberculosis
# Lineage -> from TB-profiler if species == M. tuberculosis, if no lineage detected with TB-profiler then revert to snpit
# dr_mutations -> Drug -> comma separated list of mutations
# Drug resistance summary -> based on WHO criteria
# csv columns <MDU ID>,<WGS species>,<Lineage>,<Predicted drug resistance>,<Rif (or other drug names)>,<DB version>
# json 
# {
#     MDU ID : {
#         WGS species,
#         Lineage,
#         Predicted drug resistance,
#         Rif or other drug,
#         DB version
#     }
# }

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)

def main(db, inputs, mode):
    
    isolates = []
    for i in inputs:
        tml = open_toml(i)
        isolate = list(tml.keys())[0]
        # print(list(tml.keys())[0])
        # print(tml)
        if tml[isolate]['snpit']['done'] == 'Yes' and tml[isolate]['tb_profiler']['done'] == 'Yes':
            isolates.append(isolate)
    collate_results(isolates=isolates, db=db, mode = mode)
    
db = snakemake.params.db_version
mode = snakemake.params.mode
inputs = snakemake.input

main(db = db,mode = mode, inputs = inputs)