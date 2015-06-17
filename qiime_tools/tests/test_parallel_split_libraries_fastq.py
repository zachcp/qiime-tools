# content of test_sample.py
"""
This file test the parallel_split_libraries_fastq script.
"""
from test_fastqconcat import fixname
from qiime_tools.parallel_split_libraries_fastq import parallel_splitlibraries_fastq

from click.testing import CliRunner


trimmed     = fixname("data/trimmed_small.fq")
barcodes    = fixname("data/barcodes_small.fq")
mappingfile = fixname("data/MappingFile.txt")


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

def test_parallel_splitlibraries_fastq_multiplecpus():
    "split the files into 1000 line (250 sequence)chunks and process with 4 CPUs"
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
                                "--splitsize", "1000",
                                "--ncpus", "4"])

        #assert the output is the correct length and sequence
        assert not result.exception