from __future__ import print_function
from future import standard_library
standard_library.install_aliases()

import click
from Bio import SeqIO

@click.command()
@click.option('--fastafile', type = click.File('r'), help="inputfastafile")
@click.option('--outputfasta',type = click.STRING, help="outputfile")
@click.option('--length')
def filter_by_length(fastafile, outputfasta, length):
    """
    from an input Fastafile, filter reads not matching the length parameter
    """
    with open(outputfasta, "w") as f:
        for rec in SeqIO.parse(fastafile, 'fasta'):
            if len(rec.seq) == length:
                f.write(">{}\n{}\n".format(rec.id, rec.seq))
    #records = (rec for rec in SeqIO.parse(fastafile,'fasta') if len(rec.seq) == length)
    #SeqIO.write(records, outputfasta,'fasta')
