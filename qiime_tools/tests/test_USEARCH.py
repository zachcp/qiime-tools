# content of test_sample.py
"""
This file test the parallel_split_libraries_fastq script.
"""
import os
import pandas as pd
from qiime_tools.USEARCH import UC_to_taxtable, load_ucfile

from click.testing import CliRunner

def fixname(filename):
    "alter a filename to correspond to the test directory"
    currentdir = os.path.dirname(os.path.realpath(__file__))
    return currentdir + "/" + filename


ucfile = fixname("data/test1.uc")


# Splitting Fastq Tests.
########################################################################################################################
########################################################################################################################
########################################################################################################################

def fixtargetcolumns(row):
    if row.target == "*":
        return row.query
    else:
        return row.target


def test_UC_data():
    df = load_ucfile(ucfile)
    df = df[df.rectype != "C"]
    df['target'] = df.apply(fixtargetcolumns,axis=1)
    df['sizes'] = df['query'].str.extract(';size=(\d+)')
    df['sizes'] = df['sizes'].astype(int)
    assert(df['sizes'].sum() == 24917)


def test_UC_to_taxtable():
    "run parallel_splitlibraries_fastq and check that the sequences have been concatenated."
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfile = "out.txt"
        result = runner.invoke(UC_to_taxtable,
                               ['--ucfile', ucfile,
                                '--outfile', outfile,
                                '--namehandling', "underscore2"])

        #assert the output is the correct length and sequence
        assert not result.exception
        df = pd.read_table(ucfile, header=None)
        assert(df.shape == (12,10))
        df  = pd.read_table(outfile)
        assert(df.shape == (10,2))
        samplesum = df[[1]].sum()
        assert(int(samplesum) == 24917)
