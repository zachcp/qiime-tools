"""
Test the Filter_by_length Utilits
"""

from random import choice

from Bio import SeqIO
from qiime_tools.filter_by_length import filter_by_length
from click.testing import CliRunner


def make_sequence(length=100):
    "generate a random DNA sequence of length=length"
    seq = [choice(['A','C','T','G']) for x in range(length)]
    return "".join(seq)


def test_filter_by_length():
    """
    make some fastas of variable length and filter a single one
    """
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfq = "out.fasta"
        #make a fasta file of sequences
        with open("sequences.fasta",'w') as f:
            for x in range(10):
                f.write(">sequence_{}\n{}\n".format(x, make_sequence(x)))

        result = runner.invoke(filter_by_length,
                               ['--fastafile', 'sequences.fasta',
                                '--outputfasta', outfq,
                                '--length', 5])
        assert not result.exception
        records = [rec for rec in SeqIO.parse(open(outfq,'r'),'fasta')]
        assert(len(records) == 1)


def test_filter_by_length2():
    """
    test that four of five correct-length records are returnde
    """
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfq = "out.fasta"
        #make a fasta file of sequences
        with open("sequences.fasta",'w') as f:
            for x in [5,5,5,5,6]:
                f.write(">sequence_{}\n{}\n".format(x, make_sequence(x)))

        result = runner.invoke(filter_by_length,
                               ['--fastafile', 'sequences.fasta',
                                '--outputfasta', outfq,
                                '--length', 5])

        assert not result.exception
        records = [rec for rec in SeqIO.parse(open(outfq,'r'),'fasta')]
        assert(len(records) == 4)

