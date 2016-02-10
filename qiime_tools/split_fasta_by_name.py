#from __future__ import print_function
#from future import standard_library
#standard_library.install_aliases()
#from builtins import zip
import os
from Bio import SeqIO
import click


@click.command()
@click.option('--fastafile', type=click.File('r'), prompt=True,help="name of the fasta file")
@click.option('--outdir', prompt=True, help="name of the output directory")
@click.option('--deletedir/--no-deletedir', default=False, help="whether or not to delete the directory")
def split_fasta_by_name(fastafile,outdir, deletedir):
    """
    split a fastqfile based on the names in the header

    """

    if deletedir is True:
        os.rmdir(outdir)

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    records = SeqIO.parse(fastafile,'fasta')

    for rec in records:
        sample = rec.id("_")[0]
        outlocation = "{}/{}.fasta".format(outdir, sample)
        with open(outlocation,'a') as f:
            f.write(">{} {}\n{}\n".format(rec.id, rec.description, str(rec.seq)))

    print("Splitting Complete")