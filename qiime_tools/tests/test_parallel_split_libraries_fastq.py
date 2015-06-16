# content of test_sample.py
"""
This file test the parallel_split_libraries_fastq script.
"""


from subprocess import call, check_call
from test_fastqconcat import fixname
from qiime_tools.parallel_split_libraries_fastq import parallel_splitlibraries_fastq

from click.testing import CliRunner


trimmed     = fixname("data/trimmeed.fq")
barcodes    = fixname("data/barcodes.fq")
mappingfile = fixname("data/MappingFile.txt")


# System Utility Tests
########################################################################################################################
########################################################################################################################
########################################################################################################################


def test_coreutils():
    """
    Checking that coreutils is on your system. Note: You must have coreutils >= 8.4

    Since we use the CLI of Split and require it to work as expected we require that 8.4 or greater is used.
    this will cause issues on MacOSX. This will fail if

    """
    check_call(['split', '--h'])
    check_call(['cat', '--h'])
    check_call(['zcat', '--h'])

def test_seqtk():
    """
    Checking that coreutils is on your system. Note: You must have coreutils >= 8.4

    Since we use the CLI of Split and require it to work as expected we require that 8.4 or greater is used.
    this will cause issues on MacOSX. This will fail if

    """
    assert(call('seqtk', shell=True) == 1)


# Splitting Fastq Tests.
########################################################################################################################
########################################################################################################################
########################################################################################################################


def test_parallel_splitlibraries_fastq():
    "run parallel_splitlibraries_fastq and check that the sequences have been concatenated."
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfq = "out.fq"
        result = runner.invoke(parallel_splitlibraries_fastq,
                               ["--fastq", trimmed,
                                "--barcode_fastq", barcodes,
                                "--outfile", outfq,
                                "--mappingfile", mappingfile,
                                "--barcodetype", "18",
                                "--qual_cutoff", "0",
                                "--splitsize", "100",
                                "--ncpus", "1"])

        #assert the output is the correct length and sequence
        assert not result.exception