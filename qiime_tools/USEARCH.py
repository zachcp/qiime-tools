# Functions for working with Usearch Files
#
#
#
#
import os
import click
import pandas as pd
import numpy as np
from .namehandling import processname


###################################################################
## Schemas and Data

names = ['Rec_type', 'clusternumber', 'Sequence_length', 'percent_identity',
         'strand', 'col6', 'col7', 'compressed_alignment', 'query_sequence',
         'target_sequence']

###################################################################
## Helper Functions

def openOTU(filename):
    if os.stat(filename).st_size > 0:
        df = pd.read_table(filename)
        mdf = pd.melt(df, id_vars=["#OTU ID"])
        #assert(list(mdf.columns) == ["OTU ID",  "variable", "value"])
        return mdf
    else:
        return None


def combineOTUfiles(filelist):
    dfs = (openOTU(f) for f in filelist)
    dfs = (df for df in dfs if isinstance(df, pd.DataFrame))
    big_df = pd.concat(dfs,  axis=0)
    otus = big_df.groupby(["#OTU ID","variable"]).sum().unstack().fillna(0)
    cols = [t[1] for t in otus.columns]
    otus.columns = cols
    return otus


def load_ucfile(filename):
    """

    See http://www.drive5.com/usearch/manual/opt_uc.html for
    more information about UC files.

    :param filename:
    :return: standard dataframe with usearch columns
    """
    col_names = ['rectype', 'clusternum', 'seqlength_clustsize', 'percent_ident',
                 'strand', 'nothing1','nothing2', 'compressed_algn', 'query', 'target']
    df = pd.read_table(filename, names = col_names)
    return df

def process_uc(filename, use_sizes=False):
    """
    Load a UC file and obtain the cluster information. By grouping
    on the 'query' and 'target' columns we are able to obtain the
    1-to-1 relationship between input sequences (query) and cluster
    centroids (target).

    See http://www.drive5.com/usearch/manual/opt_uc.html for
    more information about UC files.

    :param filename:
    :return: two column data frame with query and target columns representing the start and finish OTUs
    """

    def fixtargetcolumns(row):
        if row.target == "*":
            return row.query
        else:
            return row.target

    df = load_ucfile(filename)
    df['target'] = df.apply(fixtargetcolumns,axis=1)

    if use_sizes:
        df['sizes'] = df['target'].str.extract(';size=(\d+)', expand=False)
        counts = df.groupby(['query', 'target']).agg({'sizes': np.sum})
    else:
        counts = df.groupby(['query', 'target']).count().reset_index()

    return counts[['query','target']]


def returntopsums(group):
    """
    This function applies will group OTUs that belong to the
    same cluster. It returns a simple sum of the columns but replaces the index
    value with the OTUID from the OTU with the highest number of reads

    :param group:
    :return:  single row dataframe with index of max OTUID and column sums of the group
    """
    if group.shape[0] == 1:
        # return single rows
        return group.drop('target',axis=1)
    else:
        #if multiple row, return the sums

        #calculate OTU with most sums
        g1 = group.drop('target', axis=1)
        sums = g1.sum(axis=1)
        maxindex = sums.idxmax()

        #get samplesums and add the index
        samplesums = g1.sum(axis=0)
        return pd.DataFrame(samplesums, columns=[maxindex]).T

def fixindex(df):
    """ make sure integers are padded with zeros. """
    def format(x):
        x_int = int(x.replace("Seq_",""))
        return "Seq_{0:07d}".format(x_int)

    idx = list(df.index)
    df.index = map(format, idx)
    df['query'] = df.index
    return df

def load_derep(derep):
    """return a 2 column dataframe with unique reads in one columns,
       and the dereped OTU name in the next column"""

    def fix_table(x):
        if x['target_sequence'] == '*':
            return x['query_sequence']
        else:
            return x['target_sequence']
    df = pd.read_table(derep, names=names)
    df = df[ ['target_sequence','query_sequence'] ]
    df.target_sequence = df.apply(fix_table,axis=1)
    df = df.rename_axis({"query_sequence":"unique_sequence_ID","target_sequence":'derep_otu'},axis=1)
    return df[["unique_sequence_ID","derep_otu"]]

def load_97(otu97):
    """return a Dataframe of unique 97% ID Clusters without individual sequences"""
    df = pd.read_table(otu97, names=names)
    df = df[ df.Rec_type != 'N'] #eliminate outlier sequences
    df['readID'] = df.query_sequence.str.split(';').str[0]
    df['SampleCode'] = df.readID.str.split('_').str[0]
    df = df.rename_axis({"target_sequence":'otu97'},axis=1)
    return df[['SampleCode','readID','otu97'] ]

def process_secondary(ucfile,label_high, label_low):
    df = pd.read_table(ucfile, names=names)
    def add_target_sequence(row):
        """ if an OTU is not joined to a cluster then
            add a target sequence name that is the original sequence
        """
        if row['target_sequence'] == '*':
            return row['query_sequence']
        else:
            return row['target_sequence']
    df = df[df.Rec_type.isin(['S','H']) ]
    df['target'] = df.apply(add_target_sequence, axis=1)
    df = df[ ['query_sequence','target']]
    df = df.rename_axis({"target":label_low,"query_sequence":label_high}, axis=1)
    return df

###################################################################
## Public Functions

@click.command()
@click.option('--otufilename', type=click.Path(exists=True), help="name of the otufile")
@click.option('--ucfilename', type=click.Path(exists=True), help="name of of the UC file")
@click.option('--outfile', prompt=True, help="name of the new UCfile")
def merge_OTU_UCfile(otufilename, ucfilename, outfile):
    """
    This function performs a groupby operation on an OTU table, using groups (clusters)
    obtained from a Usearch OTU file. The logic here is to be able to cluster sequences using
    usearch and compress the OTU file correspondingly.

    :param otufilename: a tab-delimited otutable output, probably generated by phyloseq
    :param ucfilename: a suearch uc file
    :return: an otu file where OTUS have been clustered based on membership in the UCfile.
    """

    #load/process the OTUS and UC data
    otus = pd.read_table(otufilename)
    #handle otu files that have already been clustered
    if otus.columns[0] == 'target':
        otus = otus.set_index('target')
    #otus = fixindex(otus)
    uc_data = process_uc(ucfilename)

    # the index values of the otutable should match exactly the
    # query values of the output
    assert(len(set(otus.index).difference(set(uc_data['query']))) == 0)
    #assert(set(otus.index) == set(uc_data['query']))

    mdf = pd.merge(otus, uc_data, left_index = True, right_on="query", how="left")
    mdf = mdf.set_index("query")
    mdf = mdf.groupby("target").apply(returntopsums)
    mdf = mdf.reset_index().drop('query', axis=1).set_index('target')
    mdf.to_csv(outfile, sep="\t")


@click.command()
@click.option('--ucfile', type=click.Path(exists=True), help="name of the otufile")
@click.option('--outfile', help="name of of the UC file")
@click.option('--namehandling',
              type=click.Choice(['underscore', "underscore2","underscore3","underscore4","underscore5",
                                 'dot',"dot2","dot3", "dot4","dot5"]),
               help="how to split out the sample name from the UCfile. the options are hardcoded.")
def UC_to_taxtable(ucfile, outfile, namehandling):
    """
     the role of this script can be boiled donw to an essential coundint of
     queriy sequence that match a target.

        q1  q2   q3   q5
     t1 s1  s2   s3   s4
     t1
     t3
     t3

     where q1,q2,q3.... are the qury sequnece
           t1,t2,t2.... are the target sequences
           s1,s2,s3.... are sizes of the quer sequence.

    Calcuate the sample name and size from teh query but use the targets as rows/OTUs.

    Samplenametypes:
    1: DFD_1128.1_M03834:5:000000000-AG1GW:1:1105:222...
       returns DFD_1128.1
    2: DFD000391.r01.w0000.pAD1.M03834:4:000000000-AG1H0:1:1104:13704:13535.003791;size=9994;
        returns DFD000391.r01.w0000.

    :return:
    """

    def fixtargetcolumns(row):
        if row['target'] == "*":
            return row['query']
        else:
            return row['target']

    df = load_ucfile(ucfile)

    #remove the redundant "S" field
    df = df[df.rectype != "S"]
    df = df.reset_index()

    #fix targetnames
    df['target'] = df.apply(fixtargetcolumns,axis=1)
    #get sizes
    df['sizes'] = df['query'].str.extract(';size=(\d+)', expand=False)
    df['sizes'] = df['sizes'].astype(int)

    #getsamplenames
    df['sample'] = df['query'].apply(lambda x: processname(x, func=namehandling))

    #aggregate samples that have more than one sample per OTU
    #this can happen if more than one reads from a sample matches
    df2 = df[['target','sample','sizes']].copy()
    df3 = df2.groupby(['target','sample']).agg({"sizes": np.sum}).reset_index()
    df4 = df3.pivot('target','sample','sizes').fillna(0)
    df4.to_csv(outfile, sep="\t")


@click.command()
@click.option('--derepfile', type=click.File('r'), prompt=True,help="name of the fastq forward file")
@click.option('--otus97', type=click.File('r'), prompt=True, help="name of the 97% Uclustfile")
@click.option('--otus95', type=click.File('r'), prompt=True, help="name of the 95% Uclustfile")
@click.option('--otus90', type=click.File('r'), prompt=True, help="name of the 90% Uclustfile")
@click.option('--otus85', type=click.File('r'), prompt=True, help="name of the 85% Uclustfile")
@click.option('--otus80', type=click.File('r'), prompt=True, help="name of the 80% Uclustfile")
@click.option('--otus75', type=click.File('r'), prompt=True, help="name of the 75% Uclustfile")
@click.option('--otus55', type=click.File('r'), prompt=True, help="name of the 55% Uclustfile")
@click.option('--otus35', type=click.File('r'), prompt=True, help="name of the 35% Uclustfile")
@click.option('--otus15', type=click.File('r'), prompt=True, help="name of the 15% Uclustfile")
@click.option('--taxtable', type=click.File('w'), prompt=True, help="output the taxfile")
def UCfiles_to_taxtable(derepfile, otus97, otus95, otus90, otus85, otus80, otus75, otus55, otus35, otus15, taxtable):
    """
    This script takes a series of uclust clusterfiles and combines them to generate a master table of
    OTU membership at different percent identities. Right now the files are hardcoded so you must
    have generated these files during your cluster processing.

    Output CSV file will look like the following and will have the same number of rows as the number of non-chimeric
    sequences in your analysis (usually a few less than the number of original reads)


    unique_sequence_ID,derep_otu,SampleCode,readID,otu97,otu95,otu90,otu85,otu80,otu75,otu55,otu35,otu15
    CA06_5000213,CA06_5000213,CA06,CA06_5000213,OTU_580501,OTU_537314,OTU_416924,OTU_473588,OTU_119311,OTU_487531,OTU_415998,OTU_297299,OTU_81113
    CA06_5000511,CA06_5000213,CA06,CA06_5000213,OTU_580501,OTU_537314,OTU_416924,OTU_473588,OTU_119311,OTU_487531,OTU_415998,OTU_297299,OTU_81113
    CA06_5000562,CA06_5000213,CA06,CA06_5000213,OTU_580501,OTU_537314,OTU_416924,OTU_473588,OTU_119311,OTU_487531,OTU_415998,OTU_297299,OTU_81113
    CA06_5000589,CA06_5000213,CA06,CA06_5000213,OTU_580501,OTU_537314,OTU_416924,OTU_473588,OTU_119311,OTU_487531,OTU_415998,OTU_297299,OTU_81113
    CA06_5000599,CA06_5000213,CA06,CA06_5000213,OTU_580501,OTU_537314,OTU_416924,OTU_473588,OTU_119311,OTU_487531,OTU_415998,OTU_297299,OTU_81113
    CA06_5000681,CA06_5000213,CA06,CA06_5000213,OTU_580501,OTU_537314,OTU_416924,OTU_473588,OTU_119311,OTU_487531,OTU_415998,OTU_297299,OTU_81113
    """

    ucfiles =  [ {"file": otus95, 'label_high': 'otu97', 'label_low': 'otu95'},
                 {"file": otus90, 'label_high': 'otu95', 'label_low': 'otu90'},
                 {"file": otus85, 'label_high': 'otu90', 'label_low': 'otu85'},
                 {"file": otus80, 'label_high': 'otu85', 'label_low': 'otu80'},
                 {"file": otus75, 'label_high': 'otu80', 'label_low': 'otu75'},
                 {"file": otus55, 'label_high': 'otu75', 'label_low': 'otu55'},
                 {"file": otus35, 'label_high': 'otu55', 'label_low': 'otu35'},
                 {"file": otus15, 'label_high': 'otu35', 'label_low': 'otu15'}]

    df_derep = load_derep(derepfile)
    df = load_97(otus97)
    for uc in ucfiles:
        df2 = process_secondary( uc['file'], label_high = uc['label_high'], label_low=uc['label_low'])
        df = pd.merge(df, df2, on= uc['label_high'])

    df.to_csv(taxtable, index=False)
    return df


@click.command()
@click.option('--taxtable', type=click.File('r'), prompt=True, help="output the taxfile")
@click.option('--outfile', type=click.File('w'), prompt=True,help="name of the otu outputfile")
@click.option('--otu', type=click.STRING, default="otu95", help="name of the column to filterby")
def taxtable_to_otutable(taxtable, outfile, otu):
    """
    Takes the taxtable output from UCfiles_to_taxtable and
    produces a qiime-ready OTU table:

    #OTU ID Sample1 Sample2 ....
    otu1    0       1       ....
    otu2    0       2
    """
    df = pd.read_csv(taxtable)

    assert(otu in df.columns)
    assert("SampleCode" in df.columns)

    df = df.groupby(['SampleCode',otu]).size().reset_index() #get sample/otu stats
    df.columns = ['SampleCode',otu, 'size']
    df = df.pivot("SampleCode",otu,'size').fillna(value=0).T.reset_index() #pivot to get the correct shape
    df = df.rename_axis({otu:"#OTU ID"},axis=1) #clean
    df.to_csv(outfile,sep="\t", index=False) #write out text file
    print('Writing Taxtable in OTUTable Format')
