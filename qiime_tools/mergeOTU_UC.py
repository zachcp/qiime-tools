import pandas as pd
import numpy as np
import click

@click.command()
@click.option('--otufilename', type=click.Path(exists=True), help="name of the otufile")
@click.option('--ucfilename', type=click.Path(exists=True), help="name of of the UC file")
@click.option('--outfile', prompt=True, help="name of the new UCfile")
def merge_OTU_UCfile(otufilename, ucfilename, outfile):
    """
    This function performs a groupby operation on an OTU tale, using groups (clusters)
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


def load_ucfile(filename):
    """

    See http://www.drive5.com/usearch/manual/opt_uc.html for
    more information about UC files.

    :param filename:
    :return: standard datafrme with usearch columns
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
        df['sizes'] = df['target'].str.extract(';size=(\d+)')
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
    """

    :param df:
    :return:
    """
    def format(x):
        x_int = int(x.replace("Seq_",""))
        return "Seq_{0:07d}".format(x_int)

    idx = list(df.index)
    df.index = map(format, idx)
    df['query'] = df.index
    return df


@click.command()
@click.option('--ucfile', type=click.Path(exists=True), help="name of the otufile")
@click.option('--outfile', help="name of of the UC file")
@click.option('--samplenametype', default=1, prompt=False, help="how to split out the sample name from the UCfile. the"
																"options are hardcoded.")
def UC_to_taxtable(ucfile, outfile, samplenametype):
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

    :return:
    """

    def fixtargetcolumns(row):
        if row.target == "*":
            return row.query
        else:
            return row.target


    df = load_ucfile(filename)

    #remove the redundant "S" field
    df = df[df.rectype != "S"]

    #fix targetnames
    df['target'] = df.apply(fixtargetcolumns,axis=1)
    #get sizes
    df['sizes'] = df['query'].str.extract(';size=(\d+)')

    #getsamplenames
    if samplefunction == 1:
        #DFD_1128.1_M03834:5:000000000-AG1GW:1:1105:222...
        df['sample'] = df['query'].apply(lambda x: "_".join(x.split("_")[:2]))
    else:
        raise ValueError("Currently only a single sample naming scheme available.")

    #aggregate
    #aggregate samples that have more than one samepl per OTU
    #this can happen if more than one reads from a sample matches
    df2= df[['target','sample','sizes']].copy()
    df3 = df2.groupby(['target','sample']).agg({"sizes": np.sum})


    df3 = df.pivot('target','sample','sizes')
    df3.to_csv(outfile)

