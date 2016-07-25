import multiprocessing
from functools import partial

import click
import re
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

@click.command()
@click.option('--forward_fastq', type=click.Path(exists=True), prompt=True,help="name of the fastq forward file")
@click.option('--reverse_fastq', type=click.Path(exists=True), prompt=True, help="name of the fastq reverse file")
@click.option('--outfile', type=click.File('w'), prompt=True, help="name of the output file")
@click.option('--discard/--no-discard', default=False, help="removes paired reads where forward or reverse are shorter than keep_left or keep_right")
@click.option('--keep_left', type=click.INT, default=250, help="how much of the the forward reads should be kept.")
@click.option('--keep_right', type=click.INT, default=175, help="how much of the reverse reads should be kept")
@click.option('--ncpus',  type=click.INT, default=4, help="number of cpus to use. A little bit of parallelization \
                                                          helps. But more than a few CPUs won't get you much benefit")
@click.option('--revcomp/--no-revcomp', default=False, help="whether to reverse complement the second file")
@click.option('--spacer/--no-spacer', default=True, help="add a spacer sequence between forward and reverse")
@click.option('--spacercharacters', default="NNNNNNNNNN", help="add a spacer sequence between forward and reverse")
@click.option('--samplename', default=None, help="provide a sample anem to prefix reads with")
@click.option('--sampledelimiter', default=None, help="optional regex for processing the samplname")
def fastqconcat(forward_fastq, reverse_fastq, outfile, discard, keep_left, keep_right, ncpus, revcomp,
                spacer, spacercharacters, samplename,sampledelimiter):
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

    #for splitting the sampel from the filename
    if samplename is None:
        samplename = ""
    else:
        if sampledelimiter is not None:
            delim = sampledelimiter + ".+$"
            samplename = re.sub(delim, "", samplename)


    fastq_f = FastqGeneralIterator(open(forward_fastq,'r'))
    fastq_r = FastqGeneralIterator(open(reverse_fastq,'r'))
    fastqs = zip(fastq_f, fastq_r)

    # use partial to create a function needing only one argument
    fastqfunc = partial(process_fastq, revcomp=revcomp, discard=discard, keep_right=keep_right, keep_left=keep_left,
                        spacer=spacer, spacercharacters=spacercharacters, samplename=samplename)

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


def process_fastq(fastqs, revcomp, discard, keep_left, keep_right, spacer, spacercharacters, samplename):
    """
    take a forward and reverse fastq records, reverse complement the fastq and
    reverse the quality. Return a string on the concatenated record.
    """
    #unpack the data
    fastqf, fastqr = fastqs
    ftitle, fseq, fqual = fastqf
    rtitle, rseq, rqual = fastqr

    #if trim then check to see length of sequences are sufficient. If not return none.
    if discard:
        if len(fseq) < keep_left or len(rseq) < keep_right:
            return None

    if revcomp:
        #get reverse complement and invert the quality score
        #revcomp means the highquality data is now at the far end
        #  so slicing will occur on the left hand side counting from the right
        rseq = str(Seq(rseq).reverse_complement())
        rqual = rqual[::-1]

        if spacer:
            newseq  = fseq[:keep_left] + spacercharacters + rseq[-keep_right:]
            newqual = fqual[:keep_left] + "".join(["A" for _ in spacercharacters]) + rqual[-keep_right:]

        else:
            newseq  = fseq[:keep_left] + rseq[-keep_right:]
            newqual = fqual[:keep_left] + rqual[-keep_right:]

        return "@%s_%s\n%s\n+\n%s\n" % (samplename,ftitle, newseq, newqual)

    else:
        #put the data together
        if spacer:
            newseq  = fseq[:keep_left] + spacercharacters + rseq[:keep_right]
            newqual = fqual[:keep_left] + "".join(["A" for char in spacercharacters]) + rqual[:keep_right]

        else:
            newseq  = fseq[:keep_left] + rseq[:keep_right]
            newqual = fqual[:keep_left] + rqual[:keep_right]


        return "@%s_%s\n%s\n+\n%s\n" % (samplename, ftitle, newseq, newqual)







