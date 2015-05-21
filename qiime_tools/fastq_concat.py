from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import zip

import multiprocessing
from functools import partial

import click
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

@click.command()
@click.option('--forward_fastq', type=click.File('r'), prompt=True,help="name of the fastq forward file")
@click.option('--reverse_fastq', type=click.File('r'), prompt=True, help="name of the fastq reverse file")
@click.option('--outfile', type=click.File('w'), prompt=True, help="name of the output file")
@click.option('--keep_left', type=click.INT, default=250, help="how much of the the forward reads should be kept.")
@click.option('--keep_right', type=click.INT, default=175, help="how much of the reverse reads should be kept")
@click.option('--ncpus',  type=click.INT, default=4, help="number of cpus to use. A little bit of parallelization \
                                                          helps. But more than a few CPUs won't get you much benefit")
@click.option('--revcomp/--no-revcomp', default=False, help="whether to reverse complement the second file")
@click.option('--spacer/--no-spacer', default=True, help="add a spacer sequence between forward and reverse")
@click.option('--spacercharacters', default="NNNNNNNNNN", help="add a spacer sequence between forward and reverse")
def fastqconcat(forward_fastq, reverse_fastq, outfile, keep_left, keep_right, ncpus, revcomp,
                spacer, spacercharacters):
    """
    This script takes two fastq files and simply concatenates them to give a single 
    concatenated read as read from the forward strand. The second strand sequence is
    reverse complemented and the the quality score is reversed as in this schematic:
    
    >fasta1_F\n
    AAATTT\n
    000111\n
    
    >fasta1_R\n
    CGCATA\n
    111222\n
    
    Outputs:\n
    >fast1_F\n
    AAATTTTATGCG\n
    000111222111\n
    
    This is intended for use with MiSeq reads where the paired ends are to be used for blasting or phylogenetic comparison.
    
    """
    fastq_f = FastqGeneralIterator(forward_fastq)
    fastq_r = FastqGeneralIterator(reverse_fastq)
    fastqs = zip(fastq_f, fastq_r)

    # use partial to create a function needing only one argument
    fastqfunc = partial(process_fastq, revcomp=revcomp, keep_right=keep_right, keep_left=keep_left,
                        spacer=spacer, spacercharacters=spacercharacters)

    if ncpus==1:
        #open and process the files by shelling out each fastq pair to the pool
        results = map(fastqfunc, fastqs)
        for result in results:
            if result:
                outfile.write(result)
    else:
        # create a pool of processing nodes
        p = multiprocessing.Pool(ncpus)
        #open and process the files by shelling out each fastq pair to the pool
        results = p.imap(fastqfunc, fastqs)
        for result in results:
            if result:
                outfile.write(result)


def process_fastq(fastqs, revcomp, keep_left, keep_right, spacer, spacercharacters):
    """
    take a forward and reverse fastq records, reverse complement the fastq and
    reverse the quality. Return a string on the concatenated record.
    """
    #unpack the data
    fastqf, fastqr = fastqs
    ftitle, fseq, fqual = fastqf
    rtitle, rseq, rqual = fastqr

    if revcomp:
        #get reverse complement and invert the quality score
        #revcomp means the highquality data is now at the far end
        #  so slicing will occur on the lef thand side counting form the right
        rseq = str(Seq(rseq).reverse_complement())
        rqual = rqual[::-1]

        if spacer:
            newseq  = fseq[:keep_left] + spacercharacters + rseq[-keep_right:]
            newqual = fqual[:keep_left] + "".join(["0" for char in spacercharacters]) + rqual[-keep_right:]
            print("New Qual: {}".format(newqual))
            print(len(newseq), len(newqual))

        else:
            newseq  = fseq[:keep_left] + rseq[-keep_right:]
            newqual = fqual[:keep_left] + rqual[-keep_right:]

        return "@%s\n%s\n+\n%s\n" % (ftitle, newseq, newqual)

    else:
        #put the data together
        if spacer:
            newseq  = fseq[:keep_left] + spacercharacters + rseq[:keep_right]
            newqual = fqual[:keep_left] + "".join(["0" for char in spacercharacters]) + rqual[:keep_right]
            print("New Qual: {}".format(newqual))
            print(len(newseq), len(newqual))

        else:
            newseq  = fseq[:keep_left] + rseq[:keep_right]
            newqual = fqual[:keep_left] + rqual[:keep_right]


        return "@%s\n%s\n+\n%s\n" % (ftitle, newseq, newqual)