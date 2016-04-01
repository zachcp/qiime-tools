#from __future__ import print_function
#from future import standard_library
#standard_library.install_aliases()
#from builtins import zip
import os
from Bio import SeqIO
from .namehandling import processname
import click


@click.command()
@click.option('--fastafile', type=click.File('r'), prompt=True,help="name of the fasta file")
@click.option('--outdir', prompt=True, help="name of the output directory")
@click.option('--spliton', type=click.Choice(['underscore', "underscore2","underscore3","underscore4","underscore5",
                                              'dot',"dot2","dot3", "dot4","dot5"]))
@click.option('--deletedir/--no-deletedir', default=False, help="whether or not to delete the directory")
def split_fasta_by_name(fastafile,outdir, spliton, deletedir):
    """
    split a fastqfile based on the names in the header. splits on teh
    record ID using either an underscore or a dot character and
    taking the first member of the split list.

    :spliton: 'underscore'
        >Sample1_XXXX  -> "Sample1
    :spliton: 'dot'
        >Sample1.XXXX  -> "Sample1
    :spliton: 'underscore2
        >Sample1_Data1_XXXX  -> "Sample1_Data1
    :spliton: 'underscore4
        >Sample1_Data1_Data2_Data3_XXXX  -> "Sample1_Data1_Data2_Data3


    """

    if deletedir is True:
        os.rmdir(outdir)

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    records = SeqIO.parse(fastafile,'fasta')


    for rec in records:
        sample = processname(rec.id, func= spliton)
        outlocation = "{}/{}.fasta".format(outdir, sample)
        with open(outlocation,'a') as f:
            f.write(">{} {}\n{}\n".format(rec.id, rec.description, str(rec.seq)))

    print("Splitting Complete")