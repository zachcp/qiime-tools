#from __future__ import print_function
#from future import standard_library
#standard_library.install_aliases()
#from builtins import zip

import click
import glob
#import multiprocessing
import os
import shutil
import sys

from Bio import SeqIO
from collections import defaultdict
from subprocess import call
from cytoolz.dicttoolz import assoc
from cytoolz.functoolz import thread_first
from multiprocess import Pool

@click.command()
@click.option('--forward_fasta', prompt=True,help="name of the fasta forward file")
@click.option('--reverse_fasta', prompt=True, help="name of the fasta reverse file")
@click.option('--barcodefile', type=click.Path(exists=True), prompt=True,help="name of the barcode file")
@click.option('--barcodelength', type=click.INT, prompt=True, help="how long is the barcode")
@click.option('--outfile', prompt=True, help="output fasta file")
@click.option('--logfile', prompt=True, help="output log file")
@click.option('--max_mismatches', type=click.INT, default=1, help="maximum difference between sequence and barcode")
@click.option('--trimsize_forward',type=click.INT, default=1000)
@click.option('--trimsize_reverse',type=click.INT, default=1000)
@click.option('--includeshort/--no-includeshort', default=False)
@click.option('--spacersequence', default="NNNNNNNNNN")
@click.option('--sampleindex', type=click.INT, default=1)
@click.option('--splitsize', type=click.INT, default=100000, help="size (in lines) to split fastq")
@click.option('--ncpus', type=click.INT, default=2, help="cpus to use")
@click.option('--reverse_complement_forward/--no-reverse_complement_forward', default=False)
@click.option('--reverse_complement_reverse/--no-reverse_complement_reverse', default=True)
@click.option('--concatfirst', type=click.Choice(['forward', 'reverse']))
def demultiplex_parallel(forward_fasta, reverse_fasta, barcodefile, barcodelength, outfile,logfile, max_mismatches,
                trimsize_forward, trimsize_reverse, includeshort, spacersequence, sampleindex, splitsize,ncpus,
                         reverse_complement_forward, reverse_complement_reverse, concatfirst):
    """
    A wrapper for `demultiplex_fasta` that splitting the input fasta file and
    does the work on the pieces. Note that inputs are strings not Clikc.Files as is
    the case for the singel script.

    :param forward_fasta:
    :param reverse_fasta:
    :param barcodefile:
    :param barcodelength:
    :param outfile:
    :param logfile:
    :param max_mismatches:
    :param trimsize_forward:
    :param trimsize_reverse:
    :param includeshort:
    :param spacersequence:
    :param sampleindex:
    :return:
    """

    def makecallstring(forward_fasta,reverse_fasta, barcodefile,barcodelength,
                   outfile, logfile, max_mismatches, trimsize_forward,trimsize_reverse,
                       reverse_complement_forward, reverse_complement_reverse, concatfirst):

        if reverse_complement_forward:
            revcompstring_F = "--reverse_complement_forward"
        else:
            revcompstring_F = "--no-reverse_complement_forward"
        if reverse_complement_reverse:
            revcompstring_R = "--reverse_complement_reverse"
        else:
            revcompstring_R = "--no-reverse_complement_reverse"


        return """
        demultiplexfasta \
                --forward_fasta  {} \
                --reverse_fasta  {} \
                --barcodefile    {} \
                --barcodelength  {} \
                --outfile        {} \
                --logfile        {} \
                --max_mismatches {} \
                --trimsize_forward {} \
                --trimsize_reverse {} \
                --no-includeshort \
                {} \
                {} \
                --concatfirst={}
        """.format(forward_fasta,reverse_fasta, barcodefile,barcodelength,
                   outfile, logfile, max_mismatches, trimsize_forward,trimsize_reverse,
                   revcompstring_F,revcompstring_R, concatfirst)


    assert splitsize % 2 == 0


    splitcommand = which_split()

    #generate the split files nad check for length equivalency
    print("Splitting the Forward Fasta File, {}".format(forward_fasta))
    call("cat {} | {} -l {} - {}".format(forward_fasta, splitcommand, splitsize,"forward_"),shell=True)

    print("Splitting the Reverse Fastq File, {}".format(reverse_fasta))
    call("cat {} | {} -l {} - {}".format(reverse_fasta, splitcommand, splitsize,"reverse_"),shell=True)

    #get the names of the split files and zip them together
    split_files_forward = sorted( glob.glob("forward_*"))
    split_files_reverse = sorted( glob.glob("reverse_*"))
    split_outfiles      = sorted([x.replace("forward_","tempout_")    for x in split_files_forward])
    split_logfiles      = sorted([x.replace("forward_","logfileout_") for x in split_files_forward])

    assert len(split_files_forward) == len(split_files_reverse)

    callstrings = []
    for (f_fasta, r_fasta, out_file, log_file) in zip(split_files_forward,split_files_reverse,split_outfiles,split_logfiles):
        callstring = makecallstring(forward_fasta    = f_fasta,
                                    reverse_fasta    = r_fasta,
                                    barcodefile      = barcodefile,
                                    barcodelength    = barcodelength,
                                    outfile          = out_file,
                                    logfile          = log_file,
                                    max_mismatches   = max_mismatches,
                                    trimsize_forward = trimsize_forward,
                                    trimsize_reverse = trimsize_reverse,
                                    reverse_complement_forward = reverse_complement_forward,
                                    reverse_complement_reverse = reverse_complement_reverse,
                                    concatfirst=concatfirst)
        callstrings.append(callstring)

    #process the split files in parallel using multiprocessing
    print("Processing the Split Files in Parallel with {} cpus".format(ncpus))
    p = Pool(ncpus)
    results = p.imap_unordered(lambda x: call(x, shell=True), callstrings)
    for r in results:
        print("Processing split file.")

    print(os.listdir('.'))
    print("cleaning up the split files....")
    p.imap(os.remove, split_files_forward)
    p.imap(os.remove, split_files_reverse)


    print("Concatenating the results to {}".format(outfile))
    call("cat tempout_* > {}".format(outfile), shell=True)

    print("Concatenating the log to results to {}".format(logfile))
    call("cat logfileout_*  > {}".format(logfile), shell=True)

    print("cleaning up the temporary files....")
    p.imap(os.remove, split_outfiles)
    p.imap(os.remove, split_logfiles)

    # check output filesize is not zero which will happen
    # if something went wrong with the splitting step due to,say,
    # an error with the mapping file
    if os.path.getsize(outfile) == 0:
        os.remove(outfile)
        print("Error with Your process.. aborting...")
        raise ValueError("Outputfile of size zero indicates an issues with your qiime setup")



@click.command()
@click.option('--forward_fasta', type=click.File('r'), prompt=True,help="name of the fasta forward file")
@click.option('--reverse_fasta', type=click.File('r'), prompt=True, help="name of the fasta reverse file")
@click.option('--barcodefile', type=click.Path(exists=True), prompt=True,help="name of the barcode file")
@click.option('--barcodelength', type=click.INT, prompt=True, help="how long is the barcode")
@click.option('--outfile', type=click.File('w'), prompt=True, help="output fasta file")
@click.option('--logfile', type=click.File('w'), prompt=True, help="output log file")
@click.option('--max_mismatches', type=click.INT, default=1, help="maximum difference between sequence and barcode")
@click.option('--trimsize_forward',type=click.INT, default=1000)
@click.option('--trimsize_reverse',type=click.INT, default=1000)
@click.option('--includeshort/--no-includeshort', default=False)
@click.option('--spacersequence', default="NNNNNNNNNN")
@click.option('--sampleindex', type=click.INT, default=1)
@click.option('--includeshort/--no-includeshort', default=False)
@click.option('--reverse_complement_forward/--no-reverse_complement_forward', default=False)
@click.option('--reverse_complement_reverse/--no-reverse_complement_reverse', default=True)
@click.option('--concatfirst', type=click.Choice(['forward', 'reverse']))
def demultiplex(forward_fasta, reverse_fasta, barcodefile, barcodelength, outfile,logfile, max_mismatches,
                trimsize_forward, trimsize_reverse, includeshort, spacersequence, sampleindex,
                reverse_complement_forward, reverse_complement_reverse, concatfirst):
    """
    Demultiplexing paired Fasta files with a barcode file.

    This is intended for use with MiSeq reads where the paired ends have
    barcodes inside of the Illlumina sequencing primers. The full barcode is
    a concatenation of barcodes on  the forward and reverse read. A two-column
    barcode file must be provided that gives the barcode-to-sample relationship.

    Fasta files are expected to have been pretrimmed for quality using, for example,
    seqtk.

    This script will:

    1. check the barcode allowing for a minimal distance to the supplied barcode file.
    2. truncate forward and reverse to specified values
    3. concatenate the forward and reverse together
    4. trim by size
    5. output fasta to specified file

    """

    # get the barcode and fasta data
    barcodes   = process_barcodefile(barcodefile, barcodelength)
    fastas     = zip(SeqIO.parse(forward_fasta, 'fasta'), SeqIO.parse(reverse_fasta,'fasta'))
    fastadicts = (fasta_to_dict(fasta) for fasta in fastas)

    # get barcode information
    fastabarcodes = (check_barcode(fastadict,
                                    barcodedict=barcodes,
                                    barcodelength=barcodelength,
                                    maxdistance=max_mismatches)
                     for fastadict in fastadicts)

    #filter sizes, and reverse complement
    fastasizetruncated = (truncate_by_size(fastadict,
                                           trimsize_forward=trimsize_forward,
                                           trimsize_reverse=trimsize_reverse)
                          for fastadict in fastabarcodes)

    #iterate through and keep relevant data
    tooshortcount = 0
    badbarcodecount = 0
    errorcount = 0
    count = 0
    samplecounts = defaultdict(int)

    for result in fastasizetruncated:
        #sampledata
        #print(result)
        forward_id   = result['forward_id']
        forward_desc = result["forward_desc"]
        forward_seq  = result["forward_sequence"]
        reverse_id   = result["reverse_id"]
        reverse_desc = result["reverse_desc"]
        reverse_seq  = result["reverse_sequence"]
        sample       = result["sample"]
        barcode      = result["barcode"]
        brcd_dist    = result["barcode_distance"]
        tooshort     = result["tooshort"]

        #accounting
        count += 1
        samplecounts[sample] += 1
        if not sample: badbarcodecount += 1
        if tooshort: tooshortcount += 1

        #write sample
        def writesample(forward_seq=forward_seq, reverse_seq=reverse_seq, concatfirst=concatfirst,
                        sample=sample,forward_id=forward_id, count=count,barcode=barcode,brcd_dist=brcd_dist):

            #do reverse complement if necessary
            if reverse_complement_forward:
                forward_seq = reversecomplement(forward_seq)
            if reverse_complement_reverse:
                reverse_seq = reversecomplement(reverse_seq)

            #concat seuences in correct orientation
            if concatfirst == "forward":
                allseq = forward_seq + spacersequence + reverse_seq
            elif concatfirst == "reverse":
                allseq = reverse_seq + spacersequence + forward_seq
            else:
                raise ValueError("concatfirst must be 'forward' or 'reverse' ")

            # write out sequences
            fastaheader = "{}.{}.{:06d} barcode:{} barcodemismatches:{}".format(
                sample, forward_id, count, barcode, brcd_dist)

            outfile.write(">{}\n{}\n".format(fastaheader,allseq))

        def shouldwritesample(sample=sample,includeshort=includeshort,tooshort=tooshort,
                              brcd_dist=brcd_dist,max_mismatches=max_mismatches):
            " encapsulate sequence-writing logic in a function"

            # Only use sequences samples that have a sample
            if not sample:
                return False

            # Ignore short sequences if the flag is false
            if includeshort is False and tooshort is True:
                return False

            # Ignore sequences with barcode mismatches above the threshold
            if brcd_dist > max_mismatches:
                return False

            return True

        shouldwrite = shouldwritesample()

        if shouldwrite == True:
            writesample()

    # write out log information
    logfile.write("""
       Barcode File: {}
       Sequenced Processed: {}
       Samples Below the Length Cutoff: {}
       Samples Unassigned due to Barcodes: {}

       """.format(barcodefile, count, tooshortcount, badbarcodecount))

    for sam, cnt in samplecounts.items():
        logfile.write("Observed Counts for Sample {}: {}\n".format(sam,cnt))

    print("Finished Demultiplexing")


def reversecomplement(s):
    "reverse complement a DNA strand"
    compdict = {'A':'T','T':'A','C':'G','G':'C'}
    srev = s.upper()[::-1] #
    return "".join([compdict[char] for char in srev])

def fasta_to_dict(fasta):
    "Provide Sequence data as a dictionary"
    fastaf, fastar =  fasta
    return {"forward_id":       fastaf.id,
            "forward_desc":     fastaf.description,
            "forward_sequence": str(fastaf.seq),
            "reverse_id":       fastar.id,
            "reverse_desc":     fastar.description,
            "reverse_sequence": str(fastar.seq)}

def check_barcode(fastadict, barcodedict, barcodelength, maxdistance):
    "check for barcode and update sample data"
    samplematch = None
    hdist = 0
    halfbarcode = int(barcodelength/2)
    fseq = fastadict['forward_sequence']
    rseq = fastadict['reverse_sequence']
    barcode = fseq[:halfbarcode] + rseq[:halfbarcode]

    #check for perfect match first:
    for	sample, samplebarcode in barcodedict.items():
        if samplebarcode == barcode:
            samplematch = sample

    #if not choose closest
    if not samplematch:
        for	sample, samplebarcode in barcodedict.items():
            hdist = hamdist(samplebarcode, barcode)
            if hdist <= maxdistance:
                samplematch = sample

    # return updated values
    return thread_first(fastadict,
                        (assoc, "sample", samplematch),
                        (assoc, "barcode", barcode),
                        (assoc, "barcode_distance", hdist),
                        (assoc, "forward_sequence", fseq[halfbarcode:]),
                        (assoc, "reverse_sequence", rseq[halfbarcode:]))

def truncate_by_size(fastadict, trimsize_forward, trimsize_reverse):
    "subset sequence and indicate if short"
    fseq = fastadict['forward_sequence']
    rseq = fastadict['reverse_sequence']
    tooshort = False
    if len(fseq) < trimsize_forward:
        tooshort= True
    if len(rseq) < trimsize_reverse:
        tooshort= True

    return thread_first(fastadict,
                        (assoc, "tooshort", tooshort),
                        (assoc, "forward_sequence", fseq[:trimsize_forward]),
                        (assoc, "reverse_sequence", rseq[:trimsize_reverse]))


def process_barcodefile(file, barcodelength):
    "Take a barcode file and return the barcode"
    data = {}
    lines = open(file,'r').readlines()
    for idx, line in enumerate(lines):
        if idx > 0:
            try:
                sample, barcode  = line.split()
            except:
                sample, barcode, *othercols = line.split()

            data[sample] = barcode

    #check data
    assert(data != {})
    for k,v in data.items():
        # check barcode lengths
        if not len(v) == barcodelength:
            raise ValueError("Barcode {} is not of expected length {}".format(v, barcodelength))

    return data

def hamdist(str1, str2):
   "Count the # of differences between equal length strings str1 and str2"
   diffs = 0
   for ch1, ch2 in zip(str1, str2):
       if ch1 != ch2:
           diffs += 1
   return diffs


def which_split():
    "get aplit command using gpslit from homebrew on mac osx"

    if sys.platform == "darwin":
        #check for homebrew coreutils
        has_gsplit = shutil.which('gsplit')
        if has_gsplit is not None:
            return("gsplit")
        else:
            return("split")
    elif 'linux' in sys.platform:
        return("split")
    else:
        raise ValueError("Only Linux and MacOSX supported")

