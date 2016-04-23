# content of test_demultiplexfasta.py.py
"""
This file tests the demultiplexfasta script. demutiplexfasta two fastqfiles and outputs a fastafile where samples
have been labelled by sample id.

"""
import glob
import pytest
import os
from qiime_tools.demultiplexfasta import demultiplex, demultiplex_parallel, process_barcodefile
from Bio import SeqIO
from click.testing import CliRunner


def fixname(filename):
    "alter a filename to correspond to the test directory"
    currentdir = os.path.dirname(os.path.realpath(__file__))
    return currentdir + "/" + filename

#files
fna1 = fixname("data/fna1.fasta")
fna2 = fixname("data/fna2.fasta")
#barcodesfile = fixname("data/smallbarcodes.txt")
barcodesfile = fixname("data/smallbarcodes2.txt")
badbarcodesfile = fixname("data/smallbarcodes3.txt")

########################################################################################################################
########################################################################################################################
# Tests

def test_sequences():
    #basic fasta checking
    recs_F = [rec for rec in SeqIO.parse(open(fna1),'fasta')]
    recs_R = [rec for rec in SeqIO.parse(open(fna2),'fasta')]
    assert len(recs_F) == 5
    assert len(recs_R) == 5

def test_process_barcodefile():
    barcodes = process_barcodefile(barcodesfile, 16)
    assert(barcodes["Sample1"] == {  "barcode":         "AAAAAAAACCCCCCCC",
                                     "forward_barcode": "AAAAAAAA",
                                     "forward_spacer":  "T",
                                     "forward_primer":  "CCCTCTCTCTCCCCCCNN",
                                     "reverse_barcode": "CCCCCCCC",
                                     "reverse_spacer":  "GA",
                                     "reverse_primer":  "CCCNNNGCTCCCNNNGCT"})
    assert (barcodes["Sample2"] == {'forward_primer': 'TCNNNNTCTCTCGTGTC',
                                    'reverse_primer': 'TCTCTCTCCCTCTCCCNNNGCT',
                                    'reverse_spacer': 'GATC',
                                    'forward_barcode': 'TTTTTTTT',
                                    'reverse_barcode': 'GGGGGGGG',
                                    'barcode': 'TTTTTTTTGGGGGGGG',
                                    'forward_spacer': 'ACT'})

def test_process_barcodefile_nonuniquebarcodes():
    "badbarcodes ahs duplicate barcodes so it should trip an error"
    with pytest.raises(ValueError):
        process_barcodefile(badbarcodesfile, 16)


def test_demultiplexfasta_basic():
    "test basic demultiplexing"
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfasta   = "out.fasta"
        outlogfile = "out.log"
        result = runner.invoke(demultiplex,
                               ['--forward_fasta', fna1,
                                '--reverse_fasta', fna2,
                                '--barcodefile', barcodesfile,
                                '--barcodelength', 16,
                                '--outfile', outfasta,
                                '--logfile', outlogfile,
                                '--trimsize_forward', 5,
                                '--trimsize_reverse', 5,
                                '--concatfirst', 'forward'])

        #assert the output is the correct length and sequence
        assert result.exit_code == 0
        assert set(os.listdir(".")) == {"out.log", "out.fasta"}
        records = [r for r in SeqIO.parse(open(outfasta,'r'),'fasta')]
        assert len(records) == 2
        for rec in records:
            assert(len(rec.seq) == 20 ) # 5 + 10 spacer + 5


def test_demultiplexfasta_maxdistance():
    "change the error threshhold and observe more sequences"
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfasta   = "out.fasta"
        outlogfile = "out.log"
        result = runner.invoke(demultiplex,
                               ['--forward_fasta', fna1,
                                '--reverse_fasta', fna2,
                                '--barcodefile', barcodesfile,
                                '--barcodelength', 16,
                                '--outfile', outfasta,
                                '--logfile', outlogfile,
                                '--trimsize_forward', 5,
                                '--trimsize_reverse', 5,
                                '--max_mismatches', 3,
                                '--concatfirst', 'forward'])

        #assert the output is the correct length and sequence
        assert result.exit_code == 0
        assert set(os.listdir(".")) == {"out.log", "out.fasta"}
        records = [r for r in SeqIO.parse(open(outfasta,'r'),'fasta')]
        assert len(records) == 3
        for rec in records:
            assert(len(rec.seq) == 20 ) # 5 + 10 spacer + 5

        with runner.isolated_filesystem():
            outfasta = "out.fasta"
            outlogfile = "out.log"
            result = runner.invoke(demultiplex,
                                   ['--forward_fasta', fna1,
                                    '--reverse_fasta', fna2,
                                    '--barcodefile', barcodesfile,
                                    '--barcodelength', 16,
                                    '--outfile', outfasta,
                                    '--logfile', outlogfile,
                                    '--trimsize_forward', 5,
                                    '--trimsize_reverse', 5,
                                    '--max_mismatches', 5,
                                    '--concatfirst', 'forward'])

            # assert the output is the correct length and sequence
            assert result.exit_code == 0
            assert set(os.listdir(".")) == {"out.log", "out.fasta"}
            records = [r for r in SeqIO.parse(open(outfasta, 'r'), 'fasta')]
            assert len(records) == 4
            for rec in records:
                assert (len(rec.seq) == 20)  # 5 + 10 spacer + 5




def test_demultiplexfasta_trimsize():
    "trimming on longer sequences should eliminate one but can be retained by --includeshort"
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfasta   = "out.fasta"
        outlogfile = "out.log"
        result = runner.invoke(demultiplex,
                               ['--forward_fasta', fna1,
                                '--reverse_fasta', fna2,
                                '--barcodefile', barcodesfile,
                                '--barcodelength', 16,
                                '--outfile', outfasta,
                                '--logfile', outlogfile,
                                '--trimsize_forward', 30,
                                '--trimsize_reverse', 30,
                                '--concatfirst', 'forward'])

        #assert the output is the correct length and sequence
        assert result.exit_code == 0
        assert set(os.listdir(".")) == set(["out.log","out.fasta"])
        records = [r for r in SeqIO.parse(open(outfasta,'r'),'fasta')]
        assert len(records) == 2
        # check sequences are the correct length
        for rec in records:
            assert(len(rec.seq) == 70) # 30 + 10 spacer + 30

    with runner.isolated_filesystem():
        outfasta   = "out.fasta"
        outlogfile = "out.log"
        result = runner.invoke(demultiplex,
                               ['--forward_fasta', fna1,
                                '--reverse_fasta', fna2,
                                '--barcodefile', barcodesfile,
                                '--barcodelength', 16,
                                '--outfile', outfasta,
                                '--logfile', outlogfile,
                                '--trimsize_forward', 30,
                                '--trimsize_reverse', 30,
                                '--includeshort',
                                '--concatfirst', 'forward'])

        #assert the output is the correct length and sequence
        assert result.exit_code == 0
        assert set(os.listdir(".")) == set(["out.log","out.fasta"])
        records = [r for r in SeqIO.parse(open(outfasta,'r'),'fasta')]
        assert len(records) == 3


def test_demultiplexfasta_parallel_basic():
    "run fastqconcat and check that the sequences have been concatenated."
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfasta   = "out.fasta"
        outlogfile = "out.log"
        result = runner.invoke(demultiplex_parallel,
                               ['--forward_fasta', fna1,
                                '--reverse_fasta', fna2,
                                '--barcodefile', barcodesfile,
                                '--splitsize', 2,
                                '--barcodelength', 16,
                                '--outfile', outfasta,
                                '--logfile', outlogfile,
                                '--trimsize_forward', 5,
                                '--trimsize_reverse', 5,
                                '--no-reverse_complement_forward',
                                '--reverse_complement_reverse',
                                '--concatfirst=forward'])

        #assert the output is the correct length and sequence
        assert result.exit_code == 0
        assert set(os.listdir(".")) == {"out.log","out.fasta"}
        records = [r for r in SeqIO.parse(open(outfasta,'r'),'fasta')]
        assert len(records) == 2
        for rec in records:
            assert(len(rec.seq) == 20 ) # 5 + 10 spacer + 5

        # check that the temp files have been cleaned
        tempfiles = glob.glob("tempout_*")
        templogfiles = glob.glob("logfileout_*")
        assert(len(tempfiles) == 0)
        assert(len(templogfiles) == 0)

