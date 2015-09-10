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
@click.option('--barcodefile', type=click.Path(), prompt=True,help="name of the fastq forward file")
@click.option('--forward_fastq', type=click.File('r'), prompt=True,help="name of the fastq forward file")
@click.option('--reverse_fastq', type=click.File('r'), prompt=True, help="name of the fastq reverse file")
@click.option('--outdir', prompt=True, help="name of the output directory")
@click.option('--barcodelength', type=click.INT, prompt=True, help="how long is the barcode")
@click.option('--maxdistance', type=click.INT, prompt=True, help="maximum difference between sequence and barcode")
@click.option('--ncpus', type=click.INT, default=2, help="CPUs to use")
def demultiplex(barcodefile, forward_fastq, reverse_fastq, outdir, barcodelength, maxdistance, ncpus):
    """
    Demultiplexing paired Fastq files with split barcodes.

    This is intended for use with MiSeq reads where the paired ends have
    barcodes inside of the Illlumina sequencing primers. The full barcode is
    a concatenation of barcodes on  the forward and reverse read. A two-column
    barcode file must be provided that gives the barcode-to-sample relationship.

    """

    #lineup the fastqs
    fastq_f = FastqGeneralIterator(forward_fastq)
    fastq_r = FastqGeneralIterator(reverse_fastq)
    fastqs = zip(fastq_f, fastq_r)

    #get the barcode
    barcodes = process_barcodefile(barcodefile, barcodelength)

    #use partial to create a function needing only one argument
    fastqfunc = partial(checkbarcode, barcodes=barcodes, barcodelength=barcodelength,maxdistance=maxdistance)

    #map the partial function across fastqpairs
    p = multiprocessing.Pool(ncpus)

    results = p.map(fastqfunc, fastqs)

    for result in results:
        if result:
            sample, fastqf, fastqr = result
            assert(type(sample) == str)
            with open("{}/{}_F.fq".format(outdir,sample),'wa') as ffq:
                ffq.write(fastqf)
            with open("{}/{}_R.fq".format(outdir,sample),'wa') as rfq:
                rfq.write(fastqr)

    print("Finished Demultiplexing")

def checkbarcode(fqs, barcodes, barcodelength, maxdistance):
    "Search FQ for presence of barcode"

    #unpack the data
    fastqf, fastqr = fqs
    ftitle, fseq, fqual = fastqf
    rtitle, rseq, rqual = fastqr

    #get the barcode
    halfbarcode = barcodelength/2
    barcode = fseq[:halfbarcode] + rseq[:halfbarcode]

    #check for perfect match first:
    match = None
    for	sample, samplebarcodes in barcodes.iteritems():
        if samplebarcodes['Full'] == barcode:
            match = sample

    #if not choose closest
    if not match:
        for	sample, samplebarcodes in barcodes.iteritems():
            if hamdist(samplebarcodes['Full'], barcode) <= maxdistance:
                match = sample

    # return trimmed values
    if match:
        ffq = "@%s\n%s\n+\n%s\n" % (ftitle, fseq[halfbarcode:], fqual[halfbarcode:])
        rfq = "@%s\n%s\n+\n%s\n" % (rtitle, rseq[halfbarcode:], rqual[halfbarcode:])
        return (match, ffq, rfq)
    else:
        return None


def process_barcodefile(file, barcodelength):
    "Take a barcode file and return a nested dict of barcode info"
    data = {}
    for idx, line in enumerate(open(file,'r')):
        if idx > 0:
            sample, forward, reverse = line.split()
            data[sample] = {"Forward": forward,
                            "Reverse": reverse,
                            "Full": forward+reverse}

    #check data
    for k,v in data.iteritems():
        # check barcode lengths
        assert(len(v['Full']) == barcodelength)

    return data


def hamdist(str1, str2):
   "Count the # of differences between equal length strings str1 and str2"
   diffs = 0
   for ch1, ch2 in zip(str1, str2):
       if ch1 != ch2:
           diffs += 1
   return diffs
