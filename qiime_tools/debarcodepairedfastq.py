from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import zip


import os
import glob
import shutil

import multiprocessing
from subprocess import call
from functools import partial
from collections import defaultdict
from Bio.SeqIO.QualityIO import FastqGeneralIterator

import click

@click.command()
@click.option('--fqf', type=click.STRING, prompt=True,help="name of the fasta file")
@click.option('--fqr', type=click.STRING, prompt=True,help="name of the quality file")
@click.option('--barcodelength', type=click.INT, prompt=True,help="name of the quality file")
@click.option('--mappingfile', type=click.STRING, prompt=True, help="name of the Qiime mapping file")
@click.option('--maxmismatches', type=click.INT, default=1, help="name of the Qiime mapping file")
@click.option('--keep_unassigned', default=False, help="name of the Qiime mapping file")
@click.option('--outdir', type=click.STRING, prompt=True, help="name of the output dir")
@click.option('--splitsize', type=click.INT, default = 1000000, help="name of the output dir")
@click.option('--ncpus', type=click.INT, default=1, help="name of the output dir")
def debarcodepairedfastq(fqf, fqr, barcodelength, mappingfile, outdir, parallel,splitsize, ncpus):
    """

    A simple debarcoding script aimed at the following specific use case:


          ======================[DNA]=========================

          =========== [Forward Read]
                                        [Reverse Read]========
          ====               [Barcodes]                   ====

    Where the concatenated barcodes of the forward and reverse read are the barcode for the sample.

    This script will iterate through the F?R mate-pairs, check the barcode within a maximum distance of the reads
     and write out per-sample FastQ files to the ouput directory. This script intentionally does NO quality score
     checking as it is intended to be used upstream of the DADA2 denoising program that strictly requires the
     quality scores.
    """
    assert splitsize % 4 == 0

    # generate the split files and check for length equivalency.
    # note that the fasta and qual files are piped through seqtk to make both filetypes into
    # two-line files.
    if ncpus > 1:
        parallel = True

    if parallel:
        print("Splitting the Forward Fastq File, {}".format(fqf))
        if isgzip(fqf):
            call("zcat {} | split -l {} - {}".format(fqf, splitsize,"forward_"),shell=True)
        else:
            call("cat {} | split -l {} - {}".format(fqf, splitsize,"forward_"),shell=True)

        print("Splitting the Reverse Fastq File, {}".format(fqr))
        if isgzip(fqr):
            call("zcat {} | split -l {} - {}".format(fqr, splitsize,"reverse_"),shell=True)
        else:
            call("cat {} | split -l {} - {}".format(fqr, splitsize,"reverse_"),shell=True)

        #get the names of the split files and zip them together
        split_files_forward = sorted(glob.glob("forward_*"))
        split_files_barcode = sorted(glob.glob("barcode_*"))
        split_files_outdir = ["out_{}".format(i) for i,x in enumerate(split_files_forward)]
        number_of_files= len(split_files_forward)
        assert len(split_files_forward) == len(split_files_barcode)
        data = zip(split_files_forward, split_files_barcode, split_files_outdir, list(range(number_of_files)))

    #process the split files in parallel using multiprocessing
    if parallel:
        print("Processing the Fastq Files in Parallel with {} cpus".format(ncpus))
    else:
        print("Processing the FastQ".format(ncpus))

    p = multiprocessing.Pool(ncpus)

    handlerfunc = partial(process_split_files, mappingdict=splitlibrarycommand,
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


def process_mappingfile(mappingfile, barcodelength):
    """ assuming the mapping file conforms to Qiime's
        format guide and has passed the validate_mapping_file.py script

        return barcode: samplename dict
    """
    barcode_dict = {}
    for line in open(mappingfile, 'r'):
        linefields = line.strip().split()
        samplename = linefields[0]
        barcode = linefields[1]
        assert(len(barcode) == barcodelength)
        barcode_dict[barcode] = samplename


def process_fastqpair(fastqpair, barcodedict, barcodelength, max_mismatch, outdir):
    """

    :param fastqpair:
    :param barcodedict:
    :param barcodelength:
    :param max_mismatch:
    :return:
    """
    fqf, fqr = fastqpair
    fastq_f = FastqGeneralIterator(open(fqf, 'r'))
    fastq_r = FastqGeneralIterator(open(fqr,'r'))
    blen = barcodelength / 2
    sample = "Unassigned"

    for (f,r) in zip(fastq_f, fastq_r):
        f_name, f_seq, f_qual = f
        r_name, r_seq, r_qual = r
        barcode = f_seq[:blen] = r_seq[:blen]
        mod_barcode = None

        # check for a perfect match
        if barcode in barcodedict.keys():
            sample = barcodedict[barcode]

        # check for a single match that is within the max_value
        else:
            # make a dict of value([samples])
            # check barcodes and return a value
            hammingdists = defaultdict(list)
            for k,v in barcodedict.items():
                dist = hamdist(k, barcode)
                hammingdists[dist].append(v)
            if min(hammingdists) <= max_mismatch:
                possible_samples  = hammingdists[min(hammingdists)]
                if len(possible_samples) == 1:
                    sample = possible_samples[0]
                    barcodelist = [k for k in hammingdists if sample in hammingdists[k]]
                    if len(barcodelist) == 1:
                        mod_barcode = barcodelist[0]

        fqfout = "{}/{}_F.fastq".format(outdir, sample)
        fqrout = "{}/{}_R.fastq".format(outdir, sample)

        if mod_barcode:
            header = "{} Barcode:{} ModifiedBarcode:{}".format(sample, barcode, mod_barcode)
        else:
            header = "{} Barcode:{} ModifiedBarcode:{}".format(sample, barcode, barcode)

        with open(fqfout,'a') as f:
            f.write(">{}\n{}\n@\n{}".format(header, f_seq[blen:], f_qual[blen:]))
        with open(fqrout,'a') as f:
            f.write(">{}\n{}\n@\n{}".format(header, r_seq[blen:], r_qual[blen:]))


def hamdist(str1, str2):
    """Count the # of differences between equal length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

def isgzip(f):
    if os.path.basename(f).split(".")[-1] in ['gz', 'gzip']:
        return(True)
    else:
        return(False)

