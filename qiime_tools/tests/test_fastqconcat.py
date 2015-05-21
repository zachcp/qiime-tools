# content of test_sample.py
"""
This file test the fastqconcat script. fastconcat takes two fastqfiles and outputs an outputfiel. Our tests use
example fastq files in form the data directory and write out an output file. We then use the output file to check
whether our commands are correct.

"""
import os
from qiime_tools.fastq_concat import fastqconcat
from Bio import SeqIO
from click.testing import CliRunner


def fixname(filename):
    "alter a filename to correspond to the test directory"
    currentdir = os.path.dirname(os.path.realpath(__file__))
    return currentdir + "/" + filename

#files
fq1 = fixname("data/fq1.fq")
fq2 = fixname("data/fq2.fq")
outfile = fixname("data/out.fq")

#make sure that the output file doesn't yet exist
if os.path.exists(outfile): os.remove(outfile)

#sequences
fq1_first =  SeqIO.parse(open(fq1,'r'),'fastq').next()
fq2_first =  SeqIO.parse(open(fq2,'r'),'fastq').next()


# Tests
########################################################################################################################
########################################################################################################################
########################################################################################################################
def test_sequences():
    #basic fastq checking
    seq1 = str(fq1_first.seq)
    seq2 = str(fq2_first.seq)
    assert len(seq1) == 301
    assert len(seq2) == 301
    print fq1_first.id
    print fq2_first.id
    assert fq1_first.id == fq2_first.id

    assert seq1 == "ATAGCGACCTAGCGTACGCGATGTACACGTCCGGCTCCACAGGGACGCCCAAGGGGGTGATGAACGCGCACGGCGGCGTGGTCAACCGGCTGGCGTGGATGGAGGCCGAATACCGGCTCGGGGCCGGCGAGGCGGTGCTGCAGAAGACGCCGTACACCTTCGACGTGTCGGTGTGGGAGTTCTTCTGGCCGCTGATGACGGGCGCGCGGCTGGTGGTGGCGCGCCCCGGCGGCCACCGCGACCCCGGCTGCCTGGTCGAGACGATCGTCGCGGAGGGGATCACAACGCTTCACGTCGTCCC"
    assert seq2 == "CACTTCGAAGACGTCGCCGGTGCGGTACAGCCGCGCGCCCGGCTCGCCGAAGGGATCGGGGACGAAGCGCTCCGCCGTCTGCCCGGCGCGGTTCAGGTATCCGCGCGCCACCTGGATCCCGCCGATGTACAGCTCGCCCGCCACCCCCGGCGGCACGGGAGCCATCCTCCCATCCAGCAGGTAGGTCCGCACGTTCCCCATCGCCCGGCCGCGCGGGACGCCCCCGGATTCGCCCTCGGCGCCCCCCCCCACCGCAACCGCCACGGCGGCCTCCGCCGGGCGCGACCGGGTGGCGCGCCCC"


def test_fastqconcat_basic():
    if os.path.exists(outfile): os.remove(outfile)

    #run fastqconcat and check the output
    runner=CliRunner()
    result = runner.invoke(fastqconcat,
                           ['--forward_fastq', fq1,
                            '--reverse_fastq', fq2,
                            '--outfile', outfile,
                            '--keep_left', 301,
                            '--keep_right', 301,
                            '--ncpus', 1,
                            '--no-revcomp',
                            '--no-spacer',
                            '--spacercharacters', "NNNNNNNNNN"])

    #assert the output is the correct length and sequence
    assert not result.exception
    firstrecord = SeqIO.parse(open(outfile,'r'),'fastq').next().seq
    assert len(firstrecord) == 602
    assert str(firstrecord) == str(fq1_first.seq) + str(fq2_first.seq)
    assert str(firstrecord) == str(fq1_first.seq) + str(fq2_first.seq)


    #remove the outputfile
    os.remove(outfile)
    assert not os.path.exists(outfile)

def test_fastqconcat_revcomp():
    "Make sure the reverse complement yields the reverse complement"
    if os.path.exists(outfile): os.remove(outfile)

    #run fastqconcat and check the output
    runner=CliRunner()
    result = runner.invoke(fastqconcat,
                           ['--forward_fastq', fq1,
                            '--reverse_fastq', fq2,
                            '--outfile', outfile,
                            '--keep_left', 301,
                            '--keep_right', 301,
                            '--ncpus', 1,
                            '--revcomp',
                            '--no-spacer',
                            '--spacercharacters', "NNNNNNNNNN"])

    #assert the output is the correct length and sequence
    assert not result.exception
    firstrecord = SeqIO.parse(open(outfile,'r'),'fastq').next().seq
    assert len(firstrecord) == 602
    assert str(firstrecord) == str(fq1_first.seq) + str(fq2_first.reverse_complement().seq)

    #remove the outputfile
    os.remove(outfile)
    assert not os.path.exists(outfile)

def test_fastqconcat_trimming():
    "Keep only the first 250 characters of each"

    if os.path.exists(outfile): os.remove(outfile)

    #run fastqconcat and check the output
    runner=CliRunner()
    result = runner.invoke(fastqconcat,
                           ['--forward_fastq', fq1,
                            '--reverse_fastq', fq2,
                            '--outfile', outfile,
                            '--keep_left', 250,
                            '--keep_right', 250,
                            '--ncpus', 1,
                            '--no-revcomp',
                            '--no-spacer',
                            '--spacercharacters', "NNNNNNNNNN"])

    #assert the output is the correct length and sequence
    assert not result.exception
    firstrecord = SeqIO.parse(open(outfile,'r'),'fastq').next().seq
    assert len(firstrecord) == 500
    assert len(str(fq1_first.seq)[:250]) == 250
    assert len(str(fq2_first.seq)[:250]) == 250
    print firstrecord
    print str(fq1_first.seq)[:250]
    print str(fq2_first.seq)[:250]

    assert str(firstrecord) == str(fq1_first.seq)[:250] + str(fq2_first.seq)[:250]

    #remove the outputfile
    os.remove(outfile)
    assert not os.path.exists(outfile)


def test_fastqconcat_spacers():
    "Test that the spacer is inserted"
    if os.path.exists(outfile): os.remove(outfile)

    #run fastqconcat and check the output
    runner=CliRunner()
    spacerchars = "NNNNNNNNNN"
    result = runner.invoke(fastqconcat,
                           ['--forward_fastq', fq1,
                            '--reverse_fastq', fq2,
                            '--outfile', outfile,
                            '--keep_left', 250,
                            '--keep_right', 250,
                            '--ncpus', 1,
                            '--no-revcomp',
                            '--spacer',
                            '--spacercharacters', spacerchars])

    #assert the output is the correct length and sequence
    assert not result.exception
    firstrecord = SeqIO.parse(open(outfile,'r'),'fastq').next().seq
    assert len(firstrecord) == 500 + len(spacerchars)
    assert str(firstrecord) == str(fq1_first.seq)[:250] + spacerchars +str(fq2_first.seq)[:250]

    #remove the outputfile
    os.remove(outfile)
    assert not os.path.exists(outfile)

def test_fastqconcat_spacers_and_revcomp():
    "Check the spacer is inserted and the second sequence is reverse complemented"
    if os.path.exists(outfile): os.remove(outfile)

    #run fastqconcat and check the output
    runner=CliRunner()
    spacerchars = "NNNNNNNNNN"
    result = runner.invoke(fastqconcat,
                           ['--forward_fastq', fq1,
                            '--reverse_fastq', fq2,
                            '--outfile', outfile,
                            '--keep_left', 250,
                            '--keep_right', 250,
                            '--ncpus', 1,
                            '--revcomp',
                            '--spacer',
                            '--spacercharacters', spacerchars])

    #assert the output is the correct length and sequence
    assert not result.exception
    firstrecord = SeqIO.parse(open(outfile,'r'),'fastq').next().seq
    assert len(firstrecord) == 500 + len(spacerchars)
    assert str(firstrecord) == str(fq1_first.seq)[:250] + spacerchars +str(fq2_first.reverse_complement().seq)[-250:]

    #remove the outputfile
    os.remove(outfile)
    assert not os.path.exists(outfile)


def test_fastqconcat_basic():
    if os.path.exists(outfile): os.remove(outfile)

    #run fastqconcat and check the output
    runner=CliRunner()

    with runner.isolated_filesystem():
        outfq = "out.fq"
        result = runner.invoke(fastqconcat,
                               ['--forward_fastq', fq1,
                                '--reverse_fastq', fq2,
                                '--outfile', outfq,
                                '--keep_left', 301,
                                '--keep_right', 301,
                                '--ncpus', 1,
                                '--no-revcomp',
                                '--no-spacer',
                                '--spacercharacters', "NNNNNNNNNN"])


        #assert the output is the correct length and sequence
        assert not result.exception
        firstrecord = SeqIO.parse(open(outfq,'r'),'fastq').next().seq
        assert len(firstrecord) == 602
        assert str(firstrecord) == str(fq1_first.seq) + str(fq2_first.seq)
        assert str(firstrecord) == str(fq1_first.seq) + str(fq2_first.seq)


        #remove the outputfile
        assert not os.path.exists(outfile)