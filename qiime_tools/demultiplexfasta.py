# Code for Demultiplexing paired-end Fasta files.
#
#
#
#

import click
import csv

from Bio import SeqIO
from collections import defaultdict
from cytoolz.dicttoolz import assoc
from cytoolz.functoolz import thread_first
from schema import Schema, And, Or


###################################################################
## Schemas
barcodeSchema = Schema(
    {"barcode":         And(str, lambda s:  12 <= len(s) <= 30),
     "forward_barcode": And(str, lambda s:  6 <= len(s) <= 15),
     "forward_spacer":  And(str, lambda s:  0 <= len(s) <= 10),
     "forward_primer":  And(str, lambda s:  12 <= len(s) <= 25),
     "reverse_barcode": And(str, lambda s:  6 <= len(s) <= 15),
     "reverse_spacer":  And(str, lambda s:  0 <= len(s) <=10),
     "reverse_primer":  And(str, lambda s:  12 <= len(s) <= 25)
     })

fastadataSchema = Schema(
    {"forward_id":       str,
     "forward_desc":     str,
     "forward_sequence": str,
     "reverse_id":       str,
     "reverse_desc":     str,
     "reverse_sequence": str,
     "sample" :          Or(str, None),
     "barcode":          And(str, lambda s:  12 < len(s) < 30),
     "barcode_distance": int,
     "tooshort":         bool,
     "spacermismatch":   bool,
     })

###################################################################
## Helper Functions

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

    samplematch      = None
    barcodedata      = None
    spacermismatch   = False
    barcode_distance = 0
    halfbarcode      = int(barcodelength/2)
    fseq    = fastadict['forward_sequence']
    rseq    = fastadict['reverse_sequence']
    barcode = fseq[:halfbarcode] + rseq[:halfbarcode]

    #check for perfect match first:
    for	sample, samplebarcodedict in barcodedict.items():
        if samplebarcodedict['barcode'] == barcode:
            samplematch = sample
            barcodedata = samplebarcodedict

    #if not choose closest
    if not samplematch:
        for	sample, samplebarcodedict in barcodedict.items():
            hdist = hamdist(samplebarcodedict['barcode'], barcode)
            if hdist <= maxdistance:
                barcode_distance = hdist
                samplematch = sample
                barcodedata = samplebarcodedict

    # trim the sequences after checking the spacer sequence between the barcode and the primer
    fseq = fseq[halfbarcode:]
    rseq = fseq[halfbarcode:]

    if barcodedata is not None:
        forward_spacer = barcodedata['forward_spacer']
        reverse_spacer = barcodedata['reverse_spacer']

        if fseq.startswith(forward_spacer):
            fseq = fseq[len(forward_spacer):]
        else:
            fseq = fseq[len(forward_spacer):]
            spacermismatch = True
        if rseq.startswith(reverse_spacer):
            rseq = rseq[len(reverse_spacer):]
        else:
            rseq = rseq[len(reverse_spacer):]
            spacermismatch = True

    # return updated values
    return thread_first(fastadict,
                        (assoc, "sample", samplematch),
                        (assoc, "spacermismatch", spacermismatch),
                        (assoc, "barcode", barcode),
                        (assoc, "barcode_distance", barcode_distance),
                        (assoc, "forward_sequence", fseq),
                        (assoc, "reverse_sequence", rseq))

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
    lines = open(file,'r')
    reader = csv.reader(lines.readlines(), delimiter='\t')

    for idx, line in enumerate(reader):
        if idx > 0:
            try:
                sample, barcode, forward_barcode, forward_spacer, forward_primer, \
                reverse_barcode, reverse_spacer, reverse_primer, *othercols = line

            except:
                raise ValueError("Barcode File must have a minimum of 8 data columns")

            # validate the barcode data
            barcodedata = barcodeSchema.validate(
                            {"barcode":        barcode,
                            "forward_barcode": forward_barcode,
                            "forward_spacer":  forward_spacer,
                            "forward_primer":  forward_primer,
                            "reverse_barcode": reverse_barcode,
                            "reverse_spacer":  reverse_spacer,
                            "reverse_primer":  reverse_primer})

            data[sample] = barcodedata

    #check data
    assert(data != {})
    for k,v in data.items():
        # check barcode lengths
        if not len(v['barcode']) == barcodelength:
            raise ValueError("Barcode {}, of sample {} is not of expected length {}".format(v['barcode'], k, barcodelength))
        #check forward and reverse barcodes
        assert(v['forward_barcode'] + v['reverse_barcode'] == v['barcode'])

    barcodes = [v['barcode'] for k,v in data.items()]
    if not len(barcodes) == len(set(barcodes)):
        raise ValueError("Barcode Values are not Unique. Please Check your Barcoding File")

    return data

def hamdist(str1, str2):
   "Count the # of differences between equal length strings str1 and str2"
   diffs = 0
   for ch1, ch2 in zip(str1, str2):
       if ch1 != ch2:
           diffs += 1
   return diffs


###################################################################
## Public Functions


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

    # validate data before progressing
    fastadata = (fastadataSchema.validate(d) for d in fastasizetruncated)

    #iterate through and keep relevant data
    tooshortcount = 0
    badbarcodecount = 0
    errorcount = 0
    count = 0
    samplecounts = defaultdict(int)

    for result in fastadata:
        #sampledata
        forward_id    = result['forward_id']
        forward_desc  = result["forward_desc"]
        forward_seq   = result["forward_sequence"]
        reverse_id    = result["reverse_id"]
        reverse_desc  = result["reverse_desc"]
        reverse_seq   = result["reverse_sequence"]
        sample        = result["sample"]
        barcode       = result["barcode"]
        brcd_dist     = result["barcode_distance"]
        tooshort       = result["tooshort"]
        spacermismatch = result['spacermismatch']

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


            #concat sequences in correct orientation
            if concatfirst == "forward":
                allseq = forward_seq + spacersequence + reverse_seq
            elif concatfirst == "reverse":
                allseq = reverse_seq + spacersequence + forward_seq
            else:
                raise ValueError("concatfirst must be 'forward' or 'reverse' ")

            # write out sequences
            fastaheader = "{}.{}.{:06d} barcode:{} barcodemismatches:{} spacermismatch: {}".format(
                sample, forward_id, count, barcode, brcd_dist, str(spacermismatch))

            outfile.write(">{}\n{}\n".format(fastaheader,allseq))

        def shouldwritesample(sample=sample,includeshort=includeshort,tooshort=tooshort,
                              brcd_dist=brcd_dist,max_mismatches=max_mismatches):
            "encapsulate sequence-writing logic in a function"

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