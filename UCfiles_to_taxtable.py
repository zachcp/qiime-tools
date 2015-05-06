from __future__ import print_function
import click
import sys
import pandas as pd
import numpy as np


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


names = ['Rec_type', 'clusternumber', 'Sequence_length', 'percent_identity',
         'strand', 'col6', 'col7', 'compressed_alignment', 'query_sequence',
         'target_sequence']

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




@click.command()
@click.option('--taxtable', type=click.File('r'), prompt=True, help="output the taxfile")
@click.option('--outfile', type=click.File('w'), prompt=True,help="name of the fastq forward file")
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
    df = df.pivot("SampleCode",otu,'size').fillna(value=0).T.reset_index() #pivot to get the correct shape
    df = df.rename_axis({otu:"#OTU ID"},axis=1) #clean
    df.to_csv(outfile,sep="\t", index=False) #write out text file
    print('Writing Taxtable in OTUTable Format')
