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
@click.option('--fastq', type=click.STRING, prompt=True,help="name of the fastq file")
@click.option('--barcode_fastq', type=click.STRING, prompt=True,help="name of the barcode fastq file")
@click.option('--outfile', type=click.STRING, prompt=True, help="name of the output dir")
@click.option('--mappingfile', type=click.STRING, prompt=True, help="name of the Qiime mapping file")
@click.option('--barcodetype', type=click.STRING, prompt=True, help="filter reads below this quality score")
@click.option('--qual_cutoff', type=click.INT, default=0, help="the maximum unacceptable Phred quality score")
@click.option('--max_bad_run_length', type=click.INT, default=3, help="the most bad nts you can have before triggering truncation")
@click.option('--sequence_max_n', type=click.INT, default=10, help="the most nubmer of n's you can have")
@click.option('--splitsize', type=click.INT, default=100000, help="size (in lines) to split fastq")
@click.option('--logfile', type=click.STRING, default="split_log.txt", help="logfile name")
@click.option('--splitlibrarycommand', type=click.STRING, default="split_libraries_fastq.py", help="size (in lines) to split fastq")
@click.option('--ncpus', type=click.INT, default=4, help="number of cpus to use")
@click.option('--retain_unassigned_reads/--no-retain_unassigned_reads', type=click.BOOL, default=False, help="retain sequences which \
                        don't map to a barcode in the mapping file (sample ID will be 'Unassigned'")
@click.option('--max_barcode_errors', default=1.5, help='maximum number of errors in barcode')
def parallel_splitlibraries_fastq(fastq, barcode_fastq, outfile, mappingfile, barcodetype, qual_cutoff,
                                  max_bad_run_length, sequence_max_n, logfile,
                                  splitsize, splitlibrarycommand,ncpus,
                                  retain_unassigned_reads, max_barcode_errors):


    """
    A wrapper around Qiime's split_libraries_fastq.py

    This script will split the files, process them in parallel, and aggregate the results. Because it can run in
    parallel it us faster - MUCH faster. However this currently only supports one fastq file and one barcode file.

    Also note that the output from the split_library_log.py files are simple concatenated togehter so there is
    not great global analysis of your sequences.

    seealso: http://qiime.org/scripts/split_libraries.html
    """


    #check environment for the presence of system split and cat
    assert splitsize % 4 == 0

    try:
        check_call(['split', '--h'])
        check_call(['cat', '--h'])
        check_call(['zcat', '--h'])
    except:
        raise StandardError("coreutils not installed,out-of-date, or you are running on Windows")
    

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
                          splitsize=splitsize, max_bad_run_length=max_bad_run_length, sequence_max_n=sequence_max_n,
                          retain_unassigned_reads=retain_unassigned_reads, max_barcode_errors=max_barcode_errors)

    results = p.imap_unordered(handlerfunc, data)
    for r in results:
        print(r)

    print("cleaning up the split files....")
    p.imap(os.remove, split_files_forward)
    p.imap(os.remove, split_files_barcode)

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




def process_split_files(data,splitlibrarycommand,
                        mappingfile, qual_cutoff, barcodetype, splitsize,
                        max_bad_run_length, sequence_max_n, retain_unassigned_reads,
                        max_barcode_errors):
    """helper function for use with functional programming.
    just lets me unpack a tuple of file names"""
    fastq,barcode_fastq,outdir,number = data



    command = [splitlibrarycommand,
                   "-i", fastq,
                   "-b", barcode_fastq, 
                   "-o", outdir,
                   "-v",
                   "-m", mappingfile,
                   "-q", str(qual_cutoff),
                   "-r", str(max_bad_run_length),
                   "-n", str(sequence_max_n),
                   "--start_seq_id", str(number * (splitsize/4)),
                   "--barcode_type", barcodetype,
                   "--retain_unassigned_reads", retain_unassigned_reads,
                   "--max_barcode_errors", max_barcode_errors]


    #if discardbadwindows:
    #    command.append("--discardbadwindows")

    call(command)
    return "Finished processing a file...."

def isgzip(f):
    if os.path.basename(f).split(".")[-1] in ['gz', 'gzip']:
        return(True)
    else:
        return(False)

