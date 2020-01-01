import json,pathlib,pandas, sys



def result_path():

    result_path = pathlib.Path('results')
    if result_path:
        return result_path
    else:
        raise FileNotFoundError


def return_species(data):
    '''
    data is the dictionary from json.load
    '''
    if data["lineage"] != []:
        if "BOV" in data["lineage"][0]["lin"]:
            return "Non tuberculosis MTBC"
        else:
            return "_M. tuberculosis_"
    else:
        return "_M. tuberculosis_"

def return_lineage(data):
    '''
    data is the dictionary from json.load
    '''
    if data["lineage"] != []:
        if not "BOV" in data["lineage"][0]["lin"]:
            return f"{data['lineage'][0]['family']} ({data['lineage'][0]['lin']})"

    else:
        return f"Lineage not resolved"

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
        'Amikacin':'No mutation detected', 
        'Aminoglycosides':'No mutation detected',
        'Bedaquiline':'No mutation detected',
        'Capreomycin':'No mutation detected',
        'Ciprofloxacin':'No mutation detected', 
        'Clofazimine':'No mutation detected',
        'Cycloserine':'No mutation detected',
        'Delamanid':'No mutation detected',
        'Ethambutol':'No mutation detected', 
        'Ethionamide':'No mutation detected', 
        'Fluoroquinolones':'No mutation detected', 
        'Isoniazid':'No mutation detected',
        'Kanamycin':'No mutation detected', 
        'Levofloxacin':'No mutation detected', 
        'Linezolid':'No mutation detected', 
        'Moxifloxacin':'No mutation detected',
        'Ofloxacin':'No mutation detected', 
        'Para-aminosalicylic_acid':'No mutation detected',
        'Pyrazinamide':'No mutation detected',
        'Rifampicin':'No mutation detected'}
    
    if data['dr_variants'] != []:
        for mut in data['dr_variants']:
            drug = mut['drug']
            mutation = check_variants(mutation = mut, drug = drug)
            if mutation != '' and drug != 'streptomycin':
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
    Extensive drug-resistance predicted rif and inh AND a flouroquinolone AND one of amikacin/capreomycin/kanamycin
    '''
    # list of drigs in flq or ack for xdr
    flqlist = ['Fluoroquinolones','Levofloxacin','Moxifloxacin','Ofloxacin','Ciprofloxacin']
    acklist = ['Amikacin', 'Capreomycin', 'Kanamycin']
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
            elif d in acklist:
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
        if flq and ack: # there are mutations in a FLQ or amikacin,capreomycin,kanamycin
            resistance = 'Extensive drug-resistance predicted'
        else:
            resistance = 'Multi-drug resistance predicted'
    
    return resistance

def make_df(result_dict):

    cols_list = ['MDU ID', 'Organism identification by WGS', 'Phylogenetic lineage', 'Predicted drug resistance','Rifampicin', 'Isoniazid', 'Pyrazinamide', 'Ethambutol','Moxifloxacin','Amikacin', 'Cycloserine','Kanamycin','Capreomycin','Para-aminosalicylic_acid','Levofloxacin','Ciprofloxacin','Ofloxacin','Fluoroquinolones', 'Ethionamide','Linezolid','Clofazimine','Aminoglycosides','Bedaquiline', 'Database version used for analysis']
    
    df = pandas.DataFrame.from_dict(data = result_dict, orient = 'index')
    df = df.reset_index()
    df = df.rename(columns = {'index': 'MDU ID'})
    df = df.reindex(cols_list, axis = 'columns')

    return df

def save_results(result_dict):

    with open('tbrnr.json', 'w') as j:
        json.dump(result_dict, j)
    df = make_df(result_dict)
    df.to_csv('tbrnr.csv', index = False)

def collate_results(isolates, db):
    
    isolates = isolates.split(',')
    resultpath = result_path()
    result_dict = {}
    for isolate in isolates:
        json_path = resultpath / f"{isolate}.results.json"
        data = json.load(open(json_path))
        species = return_species(data)
        lineage = return_lineage(data)
        drugs = get_dr_variants(data)
        resistance = calculate_resistance(drugs)
        
        d = {
            'Organism identification by WGS':species,
            'Phylogenetic lineage':lineage,
            'Database version used for analysis':db,
            'Predicted drug resistance':resistance
            }
        for drug in drugs:
            d[drug] = drugs[drug]
        result_dict[isolate] = d
    
    save_results(result_dict)



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



if __name__ == "__main__":

    if len(sys.argv) == 3:
        isolates = sys.argv[1]
        print(isolates)
        db = sys.argv[2]
        collate_results(isolates=isolates, db=db)
        
