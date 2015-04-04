from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import zip


import gzip
import os
import glob
import shutil

import multiprocessing
from subprocess import call, check_call
from functools import partial

import click
from fastq_concat import isgzip


@click.command()
@click.option('--fastq', type=click.STRING, prompt=True,help="name of the fastq file")
@click.option('--barcode_fastq', type=click.STRING, prompt=True,help="name of the barcode fastq file")
@click.option('--outfile', type=click.STRING, prompt=True, help="name of the output dir")
@click.option('--mappingfile', type=click.STRING, prompt=True, help="name of the Qiime mapping file")
@click.option('--barcodetype', type=click.STRING, prompt=True, help="filter reads below this quality score")
@click.option('--qual_cutoff', type=click.INT, default=19, help="filter reads below this quality score")
@click.option('--splitsize', type=click.INT, default=100000, help="size (in lines) to split fastq")
@click.option('--splitlibrarycommand', type=click.STRING, default="split_libraries_fastq.py", help="size (in lines) to split fastq")
@click.option('--ncpus', type=click.INT, default=4, help="number of cpus to use")
def parallel_splitlibraries_fastq(fastq, barcode_fastq, outfile, mappingfile, barcodetype,qual_cutoff, splitsize, splitlibrarycommand,ncpus):
    """
    A wrapper around Qiime's split_libraries_fastq.py

    This script will split the files, process them in parallel, and aggregate the results. Because it can run in
    parallel it us faster - MUCH faster. However thi currently only supports one fastq file and one barcode file.

    Also note that the output from the split_library_log.py files are simple concatenated togehter so there is
    not great global analysis of your sequences.
    """

    #check environment for the presence of system split and cat
    assert splitsize % 4 == 0

    try:
        check_call(['split', '--h'])
        check_call(['cat', '--h'])
        check_call(['zcat', '--h'])
    except:
        raise StandardError("coreutils not installed or you are running on Windows")
    

    #generate the split files nad check for length equivalency
    print("Splitting the Fastq File, {}".format(fastq))
    if isgzip(fastq):
        call("zcat {} | split -l {} - {}".format(fastq, splitsize,"forward_"),shell=True)
    else:
        call("cat {} | split -l {} - {}".format(fastq, splitsize,"forward_"),shell=True)

    print("Splitting the Barcode Fastq File, {}".format(barcode_fastq))
    if isgzip(barcode_fastq):
        call("zcat {} | split -l {} - {}".format(barcode_fastq, splitsize,"barcode_"),shell=True)
    else:
        call("cat {} | split -l {} - {}".format(barcode_fastq, splitsize,"barcode_"),shell=True)

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
                          splitsize=splitsize)

    results = p.imap_unordered(handlerfunc, data)
    for r in results:
        print(r)

    print("cleaning up the split files....")
    p.imap(os.remove, split_files_forward)
    p.imap(os.remove, split_files_barcode)

    print("Concatenating the results to {}".format(outfile))
    call("cat out_*/seqs.fna > {}".format(outfile), shell=True)

    print("cleaning up the temporary files....")
    p.map(shutil.rmtree	, split_files_outdir)


def process_split_files(data,splitlibrarycommand, mappingfile, qual_cutoff, barcodetype, splitsize):
    """helper function for use with functional programming.
    just lets me unpack a tuple of file names"""
    fastq,barcode_fastq,outdir,number = data

    command = [splitlibrarycommand,
                   "-i", fastq,
                   "-b", barcode_fastq, 
                   "-o", outdir,
                   "-m", mappingfile,
                   "-q", str(qual_cutoff),
                   "--start_seq_id", str(number * (splitsize/4)),
                   "--barcode_type", barcodetype]

    call(command)
    return "Finished processing a file...."

