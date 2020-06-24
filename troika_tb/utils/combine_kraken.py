import toml, pathlib, subprocess, sys, pandas
from snakemake import shell

def combine_kraken(inputs):

    tab = pandas.DataFrame()

    for i in inputs:
        tml = open_toml(i)
        isolate = list(tml.keys())[0]
        if tml[isolate]['kraken']['done'] == 'Yes':
            top_3 = tml[isolate]['kraken']['top_3']
            d = {
                'Isolate':isolate,
                'Match #1': top_3[0][0],
                '%1':top_3[0][1],
                'Match #2': top_3[1][0],
                '%2':top_3[1][1],
                'Match #3': top_3[2][0],
                '%3':top_3[2][1],
            }
            df = pandas.DataFrame(d, index = [0])
            df['Quality'] = 'PASS'
        else:
            df = pandas.DataFrame(data = {'Isolate': isolate,
            'Quality': 'kmer-id not determined',
            'Match #1': '', '%1': '', 'Match #2': '', '%2': '', 'Match #3': '', '%3':''},
            index = [0])
        
        if tab.empty:
            tab = df
        else:
            tab = tab.append(df, sort = True)
    tab = tab[['Isolate', 'Match #1', '%1', 'Match #2', '%2', 'Match #3', '%3', 'Quality']]
    tab.to_csv('mtb.tab', sep= '\t', index = False)
    subprocess.run(f"sed -i 's/%[0-9]/%/g' mtb.tab", shell=True)
    return tab

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs):
    
    # set up data dict
    tab = combine_kraken(inputs)
    data = {}
    data['kraken'] = tab.to_dict(orient= 'records')
    write_toml(data = data, output= f'kraken.toml')

inputs = snakemake.input

main(inputs = inputs)