"""
This file tests the demultiplex namespace including barcode file parsing and fastQ f/r looping.

"""
import os

from qiime_tools.demultiplex import process_barcodefile, demultiplex, checkbarcode
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from click.testing import CliRunner


def fixname(filename):
    "alter a filename to correspond to the test directory"
    currentdir = os.path.dirname(os.path.realpath(__file__))
    return currentdir + "/" + filename

#files
barcodes = fixname("data/testbarcodes.txt")
fq1 = fixname("data/fq1.fq")
fq2 = fixname("data/fq2.fq")


# Tests
########################################################################################################################
########################################################################################################################
########################################################################################################################
def test_barcodefileparsing():
    #basic fastq checking
    data = process_barcodefile(barcodes, 20)
    assert(data)
    assert(type(data) ==  dict)
    assert("SampleID" not in data.keys())
    assert(data["SampleAB"]["Forward"] + data["SampleAB"]["Reverse"] == data["SampleAB"]["Full"])


def test_check_barcode():
    #check the barcode extraction the first data
    ffq = FastqGeneralIterator(open(fq1,'r'))
    rfq = FastqGeneralIterator(open(fq2,'r'))
    fqs = zip(ffq,rfq)[0]
    data = process_barcodefile(barcodes, 20)
    sample, stringf, stringr = checkbarcode(fqs, barcodes=data, barcodelength=20, maxdistance=0)
    assert(sample == "SampleAB")
    assert(stringf.startswith("@M00587:123:000000000-AEFC7:1:1109:14040:5820 1:N:0:1\nAGCGTA")) #AGCGTA is beginning after the primer
    assert(stringr.startswith("@M00587:123:000000000-AEFC7:1:1109:14040:5820 2:N:0:1\nACGTCG")) #ACGTCG is beginning after the primer

def test_demultiplex():
    #check for basic demultiplexing

    runner=CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(demultiplex, ['--barcodefile', barcodes,
                                             '--forward_fastq', fq1,
                                             '--reverse_fastq', fq2,
                                             '--outdir', ".",
                                             '--barcodelength', 20,
                                             '--maxdistance', 1,
                                             '--ncpus', 1])

        #assert the output is the correct length and sequence has the qualscore removed
        assert not result.exception
        assert os.path.exists("SampleAB_F.fq")
        assert os.path.exists("SampleAB_R.fq")
        firstrecord = SeqIO.parse(open("SampleAB_F.fq",'r'),'fastq').next().seq
        secondrecord = SeqIO.parse(open("SampleAB_R.fq",'r'),'fastq').next().seq
        assert(firstrecord == "AGCGTACGCGATGTACACGTCCGGCTCCACAGGGACGCCCAAGGGGGTGATGAACGCGCACGGCGGCGTGGTCAACCGGCTGGCGTGGATGGAGGCCGAATACCGGCTCGGGGCCGGCGAGGCGGTGCTGCAGAAGACGCCGTACACCTTCGACGTGTCGGTGTGGGAGTTCTTCTGGCCGCTGATGACGGGCGCGCGGCTGGTGGTGGCGCGCCCCGGCGGCCACCGCGACCCCGGCTGCCTGGTCGAGACGATCGTCGCGGAGGGGATCACAACGCTTCACGTCGTCCC")
        assert(secondrecord == "ACGTCGCCGGTGCGGTACAGCCGCGCGCCCGGCTCGCCGAAGGGATCGGGGACGAAGCGCTCCGCCGTCTGCCCGGCGCGGTTCAGGTATCCGCGCGCCACCTGGATCCCGCCGATGTACAGCTCGCCCGCCACCCCCGGCGGCACGGGAGCCATCCTCCCATCCAGCAGGTAGGTCCGCACGTTCCCCATCGCCCGGCCGCGCGGGACGCCCCCGGATTCGCCCTCGGCGCCCCCCCCCACCGCAACCGCCACGGCGGCCTCCGCCGGGCGCGACCGGGTGGCGCGCCCC")


def test_demultiplex_parallel():
    #check for basic demultiplexing

    runner=CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            demultiplex,
               ['--barcodefile', barcodes,
                '--forward_fastq', fq1,
                '--reverse_fastq', fq2,
                '--outdir', ".",
                '--barcodelength', 20,
                '--maxdistance', 1,
                '--ncpus', 2])

        #assert the output is the correct length and sequence
        assert result.exit_code == 0
        assert not result.exception
