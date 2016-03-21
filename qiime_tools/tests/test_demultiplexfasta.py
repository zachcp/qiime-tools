# content of test_demultiplexfasta.py.py
"""
This file tests the demultiplexfasta script. demutiplexfasta two fastqfiles and outputs a fastafile where samples
have been labelled by sample id.

"""
import os
from qiime_tools.demultiplexfasta import demultiplex, process_barcodefile
from Bio import SeqIO
from click.testing import CliRunner


def fixname(filename):
    "alter a filename to correspond to the test directory"
    currentdir = os.path.dirname(os.path.realpath(__file__))
    return currentdir + "/" + filename

#files
fna1 = fixname("data/fna1.fasta")
fna2 = fixname("data/fna2.fasta")
barcodesfile = fixname("data/smallbarcodes.txt")

########################################################################################################################
########################################################################################################################
# Tests

def test_sequences():
    #basic fasta checking
    recs_F = [rec for rec in SeqIO.parse(open(fna1),'fasta')]
    recs_R = [rec for rec in SeqIO.parse(open(fna2),'fasta')]
    assert len(recs_F) == 4
    assert len(recs_R) == 4

def test_process_barcodefile():
    barcodes = process_barcodefile(barcodesfile, 8)
    assert(barcodes["Sample1"] == "AAAAAAAA")

def test_demultiplexfasta_basic():
    "run fastqconcat and check that the sequences have been concatenated."
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfasta   = "out.fasta"
        outlogfile = "out.log"
        result = runner.invoke(demultiplex,
                               ['--forward_fasta', fna1,
                                '--reverse_fasta', fna2,
                                '--barcodefile', barcodesfile,
                                '--barcodelength', 8,
                                '--outfile', outfasta,
                                '--logfile', outlogfile,
                                '--trimsize_forward', 5,
                                '--trimsize_reverse', 5])

        #assert the output is the correct length and sequence
        assert result.exit_code == 0
        assert set(os.listdir(".")) == set(["out.log","out.fasta"])
        records = [r for r in SeqIO.parse(open(outfasta,'r'),'fasta')]
        assert len(records) == 3
        for rec in records:
            assert(len(rec.seq) == 20 ) # 5 + 10 spacer + 5


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
                                '--barcodelength', 8,
                                '--outfile', outfasta,
                                '--logfile', outlogfile,
                                '--trimsize_forward', 30,
                                '--trimsize_reverse', 30])

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
                                '--barcodelength', 8,
                                '--outfile', outfasta,
                                '--logfile', outlogfile,
                                '--trimsize_forward', 30,
                                '--trimsize_reverse', 30,
                                '--includeshort'])

        #assert the output is the correct length and sequence
        assert result.exit_code == 0
        assert set(os.listdir(".")) == set(["out.log","out.fasta"])
        records = [r for r in SeqIO.parse(open(outfasta,'r'),'fasta')]
        assert len(records) == 3

