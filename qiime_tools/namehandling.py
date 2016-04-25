# This namespace handles the names of samples that will be found
# in the fasta header.
#
#
import os
import click

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

###################################################################
## Helper Functions

def processname(name, func):
    """
    Handle the Sample name by munging the name accorging to a predefined function.
    Note that the output returns labels seperated by the 'dot' character for
    better compatability with the newick format.

    :param name: the name of a sample
    :param func: te processing function to handle the sample name
    :return: the modifiedsaple output
    """
    if func == "underscore":
        return name.split("_")[0]
    elif func == "underscore2":
        return  ".".join(name.split("_")[:2])
    elif func == "underscore3":
        return  ".".join(name.split("_")[:3])
    elif func == "underscore4":
        return  ".".join(name.split("_")[:4])
    elif func == "underscore5":
        return  ".".join(name.split("_")[:5])
    elif func == "dot":
        return  name.split(".")[0]
    elif func == "dot2":
        return  ".".join(name.split(".")[:2])
    elif func == "dot3":
        return  ".".join(name.split(".")[:3])
    elif func == "dot4":
        return  ".".join(name.split(".")[:4])
    elif func == "dot5":
        return  ".".join(name.split(".")[:5])
    else:
        raise ValueError("You are attempting to process your sample names with an undefined function")

########################################################################################
## Public Functions


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


@click.command()
@click.option('--fastq', type=click.File('r'), prompt=True,help="name of the fastq file")
@click.option('--suffix', default="", help="suffix to add to name. useful for discriminating Forward/Reverse reads")
@click.option('--outdir', prompt=True, help="name of the output directory")
def split_fastq_by_name(fastq,suffix, outdir):
    """
    split a fastqfile based on the names in the header

    """
    fastqs = FastqGeneralIterator(fastq)


    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for fq in fastqs:
        #The record name is {{name}}_{{number}}
        header, sequence, quality = fq
        sample = header.split("_")[0]
        outlocation = "{}/{}{}.fastq".format(outdir, sample, suffix)
        #print(sample, outlocation)
        with open(outlocation,'a') as f:
            f.write("@{}\n{}\n+\n{}\n".format(header,
                                            sequence,
                                            quality))
    print("Splitting Complete")