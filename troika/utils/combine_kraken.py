import pandas, pathlib, subprocess

def combine(input_list, output):
    kfiles = f"{input}".split()
    id_table = pandas.DataFrame()
    for k in kfiles:
        kraken = pathlib.Path(k)
        df = pandas.read_csv(kraken, sep = "\t", header =None, names = ['percentage', 'frag1', 'frag2','code','taxon','name'])
        df['percentage'] = df['percentage'].apply(lambda x:float(x.strip('%')) if isinstance(x, str) == True else float(x)) #remove % from columns
        df = df.sort_values(by = ['percentage'], ascending = False)
        df = df[df['code'].isin(['U','S'])]     
        df = df.reset_index(drop = True) 
        tdf = pandas.DataFrame()
        d = {'Isolate': f"{kraken.parts[0]}",    
                '#1 Match': df.ix[0,'name'].strip(), '%1': df.ix[0,'percentage'],
                '#2 Match': df.ix[1,'name'].strip(), '%2': df.ix[1,'percentage'],       
                '#3 Match': df.ix[2,'name'].strip(), '%3': df.ix[2,'percentage'] ,
                '#4 Match': df.ix[3,'name'].strip(), '%4': df.ix[3,'percentage']
                }

        tdf = pandas.DataFrame(data = d, index= [0])
        if id_table.empty:
                id_table = tdf
        else:
                id_table = id_table.append(tdf, sort = True)
    id_table = id_table[['Isolate','#1 Match' ,'%1','#2 Match','%2','#3 Match','%3']]
    id_table.to_csv(f"{output}", sep ="\t", index = False)
    subprocess.run(f"sed -i 's/%[0-9]/%/g' {output}", shell=True)



if __name__ == "__main__":

    if len(sys.argv) == 2:
        input_list = sys.argv[1]
        output= sys.argv[2]
        combine(input_list = input_list, output = output)
        