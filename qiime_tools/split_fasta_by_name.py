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
        if spliton == "underscore":
            sample = rec.id.split("_")[0]
        elif spliton == "underscore2":
            sample = ".".join(rec.id.split("_")[:2])
        elif spliton == "underscore3":
            sample = ".".join(rec.id.split("_")[:3])
        elif spliton == "underscore4":
            sample = ".".join(rec.id.split("_")[:4])
        elif spliton == "underscore5":
            sample = ".".join(rec.id.split("_")[:5])
        elif spliton == "dot":
            sample = rec.id.split(".")[0]
        elif spliton == "dot2":
            sample = ".".join(rec.id.split("_")[:2])
        elif spliton == "dot3":
            sample = ".".join(rec.id.split("_")[:3])
        elif spliton == "dot4":
            sample = ".".join(rec.id.split("_")[:4])
        elif spliton == "dot5":
            sample = ".".join(rec.id.split("_")[:5])
        else:
            raise ValueError("names are split by underscores or dots")

        outlocation = "{}/{}.fasta".format(outdir, sample)
        with open(outlocation,'a') as f:
            f.write(">{} {}\n{}\n".format(rec.id, rec.description, str(rec.seq)))

    print("Splitting Complete")