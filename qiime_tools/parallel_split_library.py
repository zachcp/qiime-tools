from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import zip


import os
import glob
import shutil

import multiprocessing
from subprocess import call, check_call
from functools import partial

import click

@click.command()
@click.option('--fasta', type=click.STRING, prompt=True,help="name of the fasta file")
@click.option('--qual', type=click.STRING, prompt=True,help="name of the quality file")
@click.option('--outfile', type=click.STRING, prompt=True, help="name of the output dir")
@click.option('--mappingfile', type=click.STRING, prompt=True, help="name of the Qiime mapping file")
@click.option('--barcodetype', type=click.STRING, prompt=True, help="filter reads below this quality score")
@click.option('--qual_cutoff', type=click.INT, default=19, help="filter reads below this quality score")
@click.option('--splitsize', type=click.INT, default=100000, help="size (in lines) to split fastq")
@click.option('--logfile', type=click.STRING, default="split_log.txt", help="logfile name")
@click.option('--splitlibrarycommand', type=click.STRING, default="split_libraries_fastq.py", help="size (in lines) to split fastq")
@click.option('--ncpus', type=click.INT, default=4, help="number of cpus to use")
@click.option('--qualwindow', type=click.INT, default=30, help="number of cpus to use")
@click.option('--barcodeerrors', type=click.INT, default=1, help="maximum allows erros in the barcode")

def parallel_split_library(fasta, qual, outfile, mappingfile, barcodetype,qual_cutoff, logfile,
                                  splitsize, splitlibrarycommand,barcodeerrors,
                                  #discardbadwindows,
                                  ncpus, qualwindow):
    """
    A wrapper around Qiime's split_libraries.py

    This script will split the files, process them in parallel, and aggregate the results. Because it can run in
    parallel it us faster - MUCH faster. However this currently only supports one fasta file and one qual file.

    Also note that the output from the split_library_log.py files are simple concatenated together so there is
    not great global analysis of your sequences.
    """
    assert splitsize % 2 == 0

    # generate the split files and check for length equivalency.
    # note that the fasta and qual files are piped through seqtk to make both filetypes into
    # two-line files.
    print("Splitting the Fasta File, {}".format(fasta))
    call("cat {}  | seqtk seq -l0 | split -l {} - {}".format(fasta, splitsize ,"forward_"),shell=True)

    print("Splitting the Qual  File, {}".format(qual))
    call("cat {}  | seqtk seq -l0 | split -l {} - {}".format(qual, splitsize,"barcode_"),shell=True)

    #get the names of the split files and zip them together
    split_files_forward = sorted(glob.glob("forward_*"))
    split_files_barcode = sorted(glob.glob("barcode_*"))
    split_files_outdir = ["out_{}".format(i) for i,x in enumerate(split_files_forward)]
    number_of_files= len(split_files_forward)
    assert len(split_files_forward) == len(split_files_barcode)
    data = zip(split_files_forward, split_files_barcode, split_files_outdir, list(range(number_of_files)))

    #process the split files in parallel using multiprocessing
    print("Processing the Split Files in Parallel with {} cpus".format(ncpus))
    p = multiprocessing.Pool(ncpus)
    handlerfunc = partial(process_split_files, splitlibrarycommand=splitlibrarycommand,
                          mappingfile=mappingfile, qual_cutoff=qual_cutoff, barcodetype=barcodetype,
                          splitsize=splitsize, qualwindow=qualwindow, barcodeerrors=barcodeerrors
                          #,discardbadwindows=discardbadwindows
                          )
    #make the output directories
    for d in split_files_outdir:
        os.mkdir(d)
    #get your results
    results = p.imap_unordered(handlerfunc, data)
    for r in results:
        print(r)

    print("cleaning up the split files....")
    p.imap(os.remove, split_files_forward)
    p.imap(os.remove, split_files_barcode)

    print("Concatenating the histogram file to results to {}".format(logfile))
    call("cat out_*/histograms.txt  > histograms.txt", shell=True)

    print("Concatenating the results to {}".format(outfile))
    call("cat out_*/seqs.fna > {}".format(outfile), shell=True)

    print("Concatenating the log to results to {}".format(logfile))
    call("cat out_*/split_library_log.txt  > {}".format(logfile), shell=True)

    print("cleaning up the temporary files....")
    p.map(shutil.rmtree	, split_files_outdir)

    # check output filesize is not zero which will happen
    # if something went wrong with the splitting step due to,say,
    # an error with the mapping file
    if os.path.getsize(outfile) == 0:
        os.remove(outfile)
        print("Error with Your process.. aborting...")
        raise ValueError("Outputfile of size zero indicates an issues with your qiime setup")


def process_split_files(data,splitlibrarycommand, mappingfile, qual_cutoff, barcodetype, splitsize,
                        qualwindow, barcodeerrors
                        #discardbadwindows
                        ):
    """helper function for use with functional programming.
    just lets me unpack a tuple of file names"""
    fasta,qual,outdir,number = data

    command = [splitlibrarycommand,
                   "-f", fasta,
                   "-q", qual,
                   "-o", outdir,
                   "-m", mappingfile,
                   "-b", str(barcodetype),
                   "-w", str(qualwindow),
                   "-e", str(barcodeerrors),
                   "-n", str(number * (splitsize/2))]
    call(command)
    return "Finished processing files {} and {}".format(fasta,qual)

def isgzip(f):
    if os.path.basename(f).split(".")[-1] in ['gz', 'gzip']:
        return(True)
    else:
        return(False)

