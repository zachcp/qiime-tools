from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import zip

import gzip
import os
import glob
import multiprocessing
from subprocess import call, check_call, Popen, PIPE
from functools import partial

import click
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

@click.command()
@click.option('--forward_fastq', type=click.STRING, prompt=True,help="name of the fastq forward file")
@click.option('--reverse_fastq', type=click.STRING, prompt=True, help="name of the fastq reverse file")
@click.option('--outfile', type=click.STRING, prompt=True, help="name of the outpuf file")
@click.option('--keep_left', type=click.INT, default=250, help="how much of the the forward reads should be kept.")
@click.option('--keep_right', type=click.INT, default=175, help="how much of the reverse reads should be kept")
@click.option('--ncpus',  type=click.INT, default=4, help="number of cpus to use. A little bit of parallelization \
                                                          helps. But more than a few CPUs won't get you much benefit")
@click.option('--gzip_out', default=True, help="whether to gzip the outputfile")
@click.option('--revcomp/--no-revcomp', default=True, help="whether to reverse complement the second file")
def fastqconcat(forward_fastq, reverse_fastq, outfile, keep_left, keep_right, ncpus, gzip_out, revcomp):
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
    concat_paired_read_files(forward_fastq=forward_fastq,reverse_fastq=reverse_fastq, 
                             outfile=outfile, keep_left=keep_left, keep_right=keep_right, 
                             ncpus=ncpus, gzip_out=gzip_out, revcomp=revcomp)

@click.command()
@click.option('--forward_fastq', type=click.STRING, prompt=True,help="name of the fastq forward file")
@click.option('--reverse_fastq', type=click.STRING, prompt=True, help="name of the fastq reverse file")
@click.option('--outfile', type=click.STRING, prompt=True, help="name of the outpuf file")
@click.option('--keep_left', type=click.INT, default=250, help="how much of the the forward reads should be kept.")
@click.option('--keep_right', type=click.INT, default=175, help="how much of the reverse reads should be kept")
@click.option('--ncpus',  type=click.INT, default=4, help="number of cpus to use. A little bit of parallelization \
                                                          helps. But more than a few CPUs won't get you much benefit")
@click.option('--splitsize',  type=click.INT, default=1000000, help="number of lines to split on.\
                                                           Must be a multiple of four.")
@click.option('--gzip_out', default=True, help="whether to gzip the outputfile")
def parallel_concat(forward_fastq, reverse_fastq, outfile, keep_left, keep_right, ncpus, splitsize, gzip_out):
    """
    This script uses the unix split command to divide the starting files and run the jobs in parallel.
    Its extra IO but the splitting and combining is done with `split` and `cat`. 
    """
    #check environment for the presence of split and cat
    assert splitsize % 4 == 0

    try:
        check_call(['split', '--h'])
        check_call(['cat', '--h'])
        check_call(['zcat', '--h'])
    except:
        raise Error("coreutils not installed or running on Windows")

    #generate the split files nad check for length equivalency
    print("Splitting the Forward Fastq File, {}".format(forward_fastq))
    if isgzip(forward_fastq):
        call("zcat {} | split -l {} - {}".format(forward_fastq, splitsize,"forward_"),shell=True)
    else:
        call("cat {} | split -l {} - {}".format(forward_fastq, splitsize,"forward_"),shell=True)

    print("Splitting the Reverse Fastq File, {}".format(reverse_fastq))
    if isgzip(reverse_fastq):
        call("zcat {} | split -l {} - {}".format(reverse_fastq, splitsize,"reverse_"),shell=True)
    else:
        call("cat {} | split -l {} - {}".format(reverse_fastq, splitsize,"reverse_"),shell=True)

    split_files_forward = sorted(glob.glob("forward_*"))
    split_files_reverse = sorted(glob.glob("reverse_*"))
    split_files_out = ["out_{}.fastq".format(i) for i,x in enumerate(split_files_forward)]
    assert len(split_files_forward) == len(split_files_reverse)
    
    #print(split_files_forward)
    #print(split_files_reverse)
    #print(split_files_out)

    triples = zip(split_files_forward, split_files_reverse, split_files_out)
    
    print("Processing the Split Files in Parallel with {} cpus".format(ncpus))
    #process the split filesin parallel
    p = multiprocessing.Pool(ncpus)
    handlerfunc = partial(process_split, keep_left=keep_left, keep_right=keep_right)
    results = p.imap(handlerfunc, triples)
    for r in results:
        print(r)
    print("cleaning up the split files....")
    p.imap(os.remove, split_files_forward)
    p.imap(os.remove, split_files_reverse)

    print("Concatenating the results to {}".format(outfile))
    call("cat out_* > {}".format(outfile), shell=True)

    print("cleaning up the temporary files....")
    p.map(os.remove, split_files_out)


def isgzip(f):
    if os.path.basename(f).split(".")[-1] in ['gz', 'gzip']:
        return(True)
    else:
        return(False)


def concat_paired_read_files(forward_fastq, reverse_fastq, outfile, keep_left, keep_right, ncpus, gzip_out, revcomp):
    """
    the meat of fastqconcat. its in a different function so it can be called internally (due to somethign funky with click)
    """
    #load the forward data
    if isgzip(forward_fastq):
        fastq_f = FastqGeneralIterator( gzip.open(forward_fastq,'r'))
    else:
        fastq_f = FastqGeneralIterator( open(forward_fastq,'r'))

    #load the reverse data
    if isgzip(reverse_fastq):
        fastq_r = FastqGeneralIterator( gzip.open(reverse_fastq,'r'))
    else:
        fastq_r = FastqGeneralIterator( open(reverse_fastq,'r'))

    #generate an outhandle
    if gzip_out:
        out_handle = gzip.open(outfile, 'w')
    else:
        out_handle = open(outfile, 'w')

    #zip the two fastq iterators together
    fastqs = zip(fastq_f, fastq_r)
    
    # use partial to create a function needing only one argument
    fastqfunc = partial(process_fastq, revcomp=revcomp, keep_right=keep_right, keep_left=keep_left)
    
    if ncpus==1:
        #open and process the files by shelling out each fastq pair to the pool
        with open(outfile, "w") as out_handle:
            results = map(fastqfunc, fastqs)
            for result in results:
                if result:
                    out_handle.write(str(result))
    else:            
        # create a pool of processing nodes
        p = multiprocessing.Pool(ncpus)   
        #open and process the files by shelling out each fastq pair to the pool
        with open(outfile, "w") as out_handle:
            results = p.imap(fastqfunc, fastqs)
            for result in results:
                if result:
                    out_handle.write(str(result))


def process_fastq(fastqs, revcomp ,keep_left, keep_right):
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
        rseq = str(Seq(rseq).reverse_complement())
        rqual = rqual[::-1]

    #put the data together
    newseq  = fseq[:keep_left] + rseq[-keep_right:]
    newqual = fqual[:keep_left] + rqual[-keep_right:]

    return "@%s\n%s\n+\n%s\n" % (ftitle, newseq, newqual)

def process_split(triple, keep_left,keep_right):
    """helper function for use with functional programming.
    just lets me unpack a tuple of file names"""
    f,r,out = triple
    concat_paired_read_files(forward_fastq=f, reverse_fastq=r, outfile=out,
                keep_left=keep_left,keep_right=keep_right, ncpus=1, gzip_out=False)
    return "Finished processing {}".format(triple)
    
