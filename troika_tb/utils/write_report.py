
import svgwrite, re, sys, subprocess, toml
from Bio import Phylo
import jinja2, pathlib, pandas, numpy, re
from packaging import version
import datetime

from snakemake import shell

class Report:
    '''
    A class to generate the tables and figures for use in the report.html
    '''
    def write_tables(self,table):
        '''
        Write a table, given a tab delimited file generate a html string
        '''
        # TODO add class isolate id to <tr>
        # TODO add class distances-isolate to tr if matrix and head-isolate to head td
        path = pathlib.Path(f"{table}")
        if path.exists():
            data = open(path).readlines()
            # for header
            header = data[0].split('\t') 
            if 'distances.tab' in table:
                tablehead = [f"<th class='{column}-head'>{column}</th>" for column in header]
            else:
                tablehead = [f"<th>{column}</th>" for column in header]
            # for body seqtablebody
            body = []
            for i in range(1,len(data)):
                raw = data[i].split('\t')  
                if 'summary_table.tab' in table:
                    row = [f"<tr class='{raw[0]} tiplab'>"]
                elif 'distances.tab' in table:
                    row = [f"<tr class='distances-{raw[0]}'>"]
                elif 'mtb' in table:
                    row = [f"<tr class='{raw[0]}-species-identification'>"]
                elif 'core_genome.tab' in table:
                    row = [f"<tr class='{raw[0]}-core-genome'>"]
                elif 'troika.tab' in table:
                    row = [f"<tr class='{raw[0]}-drug-resistance'>"]
                elif 'seqdata.tab' in table:
                    row = [f"<tr class='{raw[0]}-sequence-data'>"]
                else:
                    row = [f"<tr>"] # TODO add class isolate id to <tr>
                if 'distances.tab' in table:
                    for d in range(len(raw)):
                        row.append(f"<td align=\"center\" class = \"{raw[0]}_{header[d]}\">{raw[d]}</td>")    
                else:
                    for d in raw:
                        row.append(f"<td align=\"center\">{d}</td>")
                row.append(f"</tr>")
                body = body + row
            return('\n'.join(tablehead),'\n'.join(body))

    

    def get_table_data(self,td):
        '''
        input:
            :reportdir: the directory where report files are kept
            :td: the dictionary for data
        output:
            :td: an updated dictionary
        '''

        for tabletype in range(len(td)):
            table = td[tabletype]['file']
            td[tabletype]['head'], td[tabletype]['body'] = self.write_tables(table=table)
        
        return(td)

    def plot_histogram(self,series, xlabel, ylabel, bins, color = "#3973ac", width = 1000, plot_height = 200):
        '''
        generic method for generating a histogram from a dataframe using bokeh
        input:
            :series:the column of the dataframe to be used
            :xlabel: label for x-axis
            :ylabel: label for y-axis
            :bins: the size of the histogram bins
            :color: color of the bars
            :width: width of the graph
            :height: height of the graph
        output:
            :script: the javascript string for insert into html
            :div: the html div for inster in html doc
        '''
        # generate histogram
        hist, edges = numpy.histogram(series, density=True, bins=bins)
        # generate bokeh plot
        p = figure(plot_width=width, plot_height=plot_height)
        p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color=color, color = color)
        # style
        p.xaxis.axis_label = xlabel
        p.yaxis.axis_label = ylabel
        p.sizing_mode = "scale_width"
        # save as png
        # export_png(p, filename="snpdensity.png")
        # get div and script and and add to dict
        script, div = components(p)
        # show(p)
        return(script,div)

    def plot_snpdensity(self,idx_file):

        '''
        generate a snp-density accross the genome plot - using the core.tab file
        input:
            :reportdir: the directory where files are kept
        out:
            :snpdensityscript: the javascript string for insert into html
            :snpdenstiydiv: the html div for inster in html doc
        '''
        # helper functions for getting the data into the right format.
        def adjust_offset(row, d):
            if row['CHR'] in d:
                return(int(d[row['CHR']]) + int(row['POS']))
            else:
                return(int(row['POS']))

        def generate_dict(idx_file):
            d = {}
            pos = 0
            offset = 0
            with open(f"{idx_file}") as f:
                for line in f:
                    l = line.split()
                    d[l[0]] = offset
                    offset = offset + int(l[1])
            return(d)
        
        # open fai file and generate the dictionary
        
        d = generate_dict(idx_file)
        # get the core file
        core = 'core.tab'
        df = pandas.read_csv(core, sep = '\t')
        # get a list of isolate names
        names = list(df.columns[3:len(df.columns)])
        # if the there is no snp in the isolate (ie same as ref then mak na - then easy to drop)
        for i in names:
                df[i]=numpy.where(df['REF'] == df[i], numpy.nan, df[i])
        # generate the offset value
        df['POS_OFFSET'] = df[['CHR', 'POS']].apply(lambda x:adjust_offset(x,d), axis = 1)
        # collect positions to get allow for histogram and dropna (no snp)
        melted_df = pandas.melt(df, id_vars=['POS_OFFSET'], value_vars=names)
        melted_df = melted_df.dropna()
        melted_df = melted_df.sort_values(by= ['POS_OFFSET'])
        # generate histogram
        # snpdensityscript, snpdensitydiv = self.plot_histogram(series=melted_df['POS']/1000,xlabel="Genome Position (MB)", ylabel="SNPS",bins=10000)
        contig_breaks = [d[value] for value in d]
        
        # return dictionary
        return(list(melted_df['POS_OFFSET']/1000))

    def plot_distances(self):

        '''
        generate a snp-density plot - using the distacnes.tab file
        input:
            :reportdir: the directory where files are kept
        out:
            :distancescript: the javascript string for insert into html
            :distancesdiv: the html div for inster in html doc
        '''
        distance = 'distances.tab'

        df = pandas.read_csv(distance, sep = '\t')
        # get a list of isolate names
        names = list(df.columns[1:len(df.columns)])
        col1 = df.columns[0]
        
        # if the there is no snp in the isolate (ie same as ref then mak na - then easy to drop)
        
        # collect positions to get allow for histogram and dropna (no snp)
        melted_df = pandas.melt(df, id_vars=[col1], value_vars=names)
        melted_df = melted_df[melted_df[col1]!= melted_df['variable']]
        # generate histogram
        # distancescript, distancediv = self.plot_histogram(series=melted_df['value'], xlabel="SNP distances", ylabel="Frequency", bins=100)
        # td['pairwisedistance'] = {'script':distancecript, 'div':distancediv}
        # return dictionary
        return(list(melted_df['value']))


    def get_tree_string(self):
        '''
        Generate a tree image from a newick
        input:
            :reportdir: the directory where report files are stored
        output:
            string reporesentation of the path to the tree image
        '''
        # get tree

        with open(f"core.treefile", 'r') as t:
            tree = t.read().strip()

        return tree

    def get_software_versions(self, software):

        '''
        Given the name of the software, find the version
        input:
            :software: the name of the software
        output:
            a string in the form of 'Name_of_Sofware v.X.Y.Z'
        '''

        version_pat = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+))?\b')

        if software == 'snp-dists':
            vs = '-v'
        else:
            vs = '--version'
        cmd = f"{software} {vs} 2>&1"
        
        p = subprocess.run(cmd, shell = True, capture_output=True, encoding = "utf-8")
        sft = p.stdout
        
        v = version_pat.search(sft)
        v = v.group()
        sft_version = f"{software} v.{v}"
        return(sft_version)
    
    def make_dict_versions(self, tools):
        '''
        Called by get_software_file to make a dictionary of tools used
        input:
            :tools: a list of tools
        output:
            a dictionary with tools and versions.
        '''
        tool_dict = {}
        for t in tools:
               v = self.get_software_versions(t)
               tool_dict[t] = v
        return(tool_dict)

    def get_software_file(self):
        '''
        get the versions of software on the system at completion of the piepline
        input:
            :reportdir: the directory where report files are stored
            :pipeline: the type of pipeline
            :assembler: the assembler used in the pipeline
        '''
        tools = ['snippy', 'snippy-core', 'snp-dists', 'iqtree', 'tb-profiler']
        
        tool_dict = self.make_dict_versions(tools)
        
        versions = ['Software versions']
        for t in tool_dict:
            versions.append(tool_dict[t])
        
        p = pathlib.Path('software_versions.tab')

        p.write_text('\n'.join(versions))

    def merge_dfs(self,start, added):
        
        if start.empty:
            start = added
        else:
            start = start.merge(added)
        return(start)

    def generate_summary(self):
        '''
        function to generate a summary table
        '''
        tabs = ['troika.tab', 'mtb.tab', 'distances.tab', 'seqdata.tab','core_genome.tab','core.treefile','core.tab']
        summary_df = pandas.DataFrame()
        df_list = []
        print(tabs)
        for tab in tabs:
            # print(df)
            if tab == 'mtb.tab' and pathlib.Path('mtb.tab').exists():
                species = pandas.read_csv(tab, sep = '\t')
                print(species)
                species = species[['Isolate', 'Match #1']]
                summary_df = self.merge_dfs(summary_df, species)
                summary_df = summary_df.rename(columns={'#1 Match': 'Species'})
            elif 'seqdata' in tab:
                seq = pandas.read_csv(pathlib.Path(tab), sep = '\t')
                seq = seq[['Isolate', 'Estimated depth']]
                summary_df = self.merge_dfs(summary_df, seq)
            elif 'core_genome' in tab:
                core = pandas.read_csv(pathlib.Path(tab), sep = '\t')
                core = core[['Isolate', '% USED']]
                summary_df = self.merge_dfs(summary_df,core)
            elif 'troika' in tab:
                res = pandas.read_csv(pathlib.Path(tab), sep = '\t')
                # print(res.columns)
                res = res[['Isolate', 'Predicted drug resistance', 'Phylogenetic lineage']]
                # res = res.rename(columns = {'MDU ID': 'Isolate'})
                summary_df = self.merge_dfs(summary_df,res)
        if 'Species' not in list(summary_df.columns):
            summary_df['Species'] = 'Species detection not performed'
        summary_file = 'summary_table.tab'
        summary_df.to_csv(summary_file, sep = '\t', index = False)

        return len(summary_df['Isolate'])

    def get_metadata(self):

        resdata = pandas.read_csv('troika.tab', sep = '\t') 
        data = {}
        for row in resdata.iterrows():
            print(row[1])
            rif =  "#adf0ad" if row[1]["Rifampicin"]== "No mutation detected" else "#f03329"
            inh = "#adf0ad" if row[1]["Isoniazid"]== "No mutation detected" else "#f03329"
            pza = "#adf0ad" if row[1]["Pyrazinamide"]== "No mutation detected" else "#f03329"
            emb = "#adf0ad" if row[1]["Ethambutol"]== "No mutation detected" else "#f03329"
            data[row[1]["Isolate"]] = {
                'Rif':{'colour':rif},
                'Inh':{'colour':inh},
                'Pza':{'colour':pza},
                'Emb':{'colour':emb}
            }
        return data,len(list(resdata['Isolate']))

    def open_toml(self,tml):

        data = toml.load(tml)

        return data
    
    def write_toml(self,data, output):
        
        with open(output, 'wt') as f:
            toml.dump(data, f)

    def main(self,workdir, resources, amr_only, idx):
        '''
        main function of the report class ties it all together
        input:
            :workdir: job directory
            :resources: the directory where templates are stored
            :amr_only: whether or not amr_only is being run
        '''
        
        
        # set up paths variables
        p = pathlib.Path('.')
        date = datetime.datetime.today().strftime("%d/%m/%y")
        # path to report data
        reportdir = pathlib.Path(workdir)
        reporthtml = pathlib.Path('index.html')
        # path to html template
        indexhtml = pathlib.Path(resources,'index.html') # replace with template 
        # newick string
        newick_path = 'core.treefile'
        newick_string = open(newick_path).read().strip()
        # save tool table
        self.get_software_file()
        
        # table dictionary for each option
        
        td = [{'file':'seqdata.tab', 'title':'Sequence Data', 'link': 'sequence-data', 'type' : 'table'}, {'file':'summary_table.tab','title':'Summary', 'link':'summary', 'type':'summary'}, {'file':'troika.tab', 'title':'Drug Resistance', 'type':'table', 'link':'drug-resistance'}]
        sp = 'mtb.tab'
        
        # TODO edit links to be title lower case separated by a -
        core_genome_td = {'file': 'core_genome.tab', 'title': 'Core Genome', 'link':'core-genome', 'type':'table'}
        snp_density_td = {'title': 'SNP density', 'link':'snp-density', 'type':'graph'}
        core_phylogeny_td = {'title': 'Phylogeny', 'link':'phylogeny', 'file': 'core.treefile', 'type': 'tree'}
        snp_distance_td = {'file': 'distances.tab', 'title':'SNP distances', 'type':'matrix', 'link':'snp-distances'}
        # list of snp tasks
        s_td = [core_genome_td,snp_density_td,core_phylogeny_td, snp_distance_td]
        # species_id_td = 
        
        # list of assembly tasks       
    
        td.extend(s_td)
        tables =['drug-resistance', 'core-genome', 'snp-distances', 'sequence-data']
        modaltables =['drug-resistance','sequence-data', 'core-genome']
        display = f"display:inline;"
    
        tables.append('versions')
        # get versions of software
        versions_td = {'file': 'software_versions.tab', 'title': 'Tools', 'type': 'versions', 'link':'versions'}
        td.append(versions_td)
        snpdistances = ''
        snpdensity = ''
        # add data to sections
        # print(td)
        for t in range(len(td)):
            # print(t)
            # TODO if table add a modal modal + link and link will be title lowercase with hyphen
            if td[t]['type'] == 'table':
                td[t]['head'], td[t]['body'] = self.write_tables(table=td[t]['file'])
            if td[t]['type'] == 'tree':
                td[t]['image'] = self.get_tree_string()
            if td[t]['link'] == 'snp-distances':
                td[t]['head'], td[t]['body'] = self.write_tables(table=td[t]['file'])
                snpdistances = self.plot_distances()
            if td[t]['link'] == 'snp-density':
                snpdensity= self.plot_snpdensity(idx_file = idx)
            if td[t]['type'] == 'versions':
                td[t]['head'], td[t]['body'] = self.write_tables(table=td[t]['file'])
            if td[t]['type'] == 'summary':
                self.generate_summary()
                td[t]['head'], td[t]['body'] = self.write_tables(table = td[t]['file'])
        metadata_string,height  = self.get_metadata()
        data = {'newick': newick_string, 
                'display': display,
                'tables': tables,
                'td': td, 
                'snpdistances':snpdistances, 
                'snpdensity': snpdensity, 
                'modaltables' : modaltables, 
                'date': date, 
                'metadata_string': metadata_string, 
                'height': height*50}

        self.write_toml(data, 'report.toml')
        # fill template
        
        report_template = jinja2.Template(pathlib.Path(indexhtml).read_text())
        reporthtml.write_text(report_template.render(data))

resources = snakemake.params.template_path
workdir = snakemake.params.wkdir
amr_only = snakemake.params.amr_only
idx = snakemake.params.idx

R = Report()
R.main(resources = resources, workdir = workdir, amr_only = amr_only, idx = idx)