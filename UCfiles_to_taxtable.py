from __future__ import print_function
import click

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
def UCfiles_to_taxtable(derepfile, otus97, otus95, otus90, otus85, otus80, otus75, otus55, otus35, otus15, outfile):
    """
    This script takes a series of uclust clusterfiles and cobmines them to generate a master table of
    OTU membership at different percent identities. Right now the files are hardcoded so you must
    have generated these files during your cluster processing.
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
    df = pd.merge(df_derep,df, left_on='derep_otu', right_on='readID')

    for uc in ucfiles:
        df2 = process_secondary( uc['file'], label_high = uc['label_high'], label_low=uc['label_low'])
        df = pd.merge(df, df2, on= uc['label_high'])

    df.to_csv(taxtable)
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
        """ if an OTU is not joined to a cluter then
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
