# content of test_sample.py
"""
This file test the parallel_split_libraries_fastq script.
"""
import pandas as pd
from test_fastqconcat import fixname
from qiime_tools.mergeOTU_UC import UC_to_taxtable

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
                                '--outfile', outfile])

        #assert the output is the correct length and sequence
        assert not result.exception
        df = pd.read_table(ucfile, header=None)
        assert(df.shape == (12,10))
        df  = pd.read_table(outfile)
        assert(df.shape == (10,2))
        assert(df[[1]].sum() == 24917)
