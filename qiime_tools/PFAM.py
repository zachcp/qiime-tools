# This Namespace covers basic utilities tha interact with the PFAM dataset.
# Most notably, this includes obtaining the nucleic acid sequences belonging
# to the proteins in a PFAM alignment By getting info as follows
#  PFAM Domain -> Uniprot Protein IDs -> ENA Nucleic Acid IDs -> ENA
#
import os
import requests
import pandas as pd
import uniprot as uni
from collections import defaultdict
from Bio import SeqIO

def getUniprotIDsfromPFAM(PFAMID, keepfiles=True):
    """
    Take a PFAM ID and obtain the Nucleotide Sequences of the PFAMs in that alignment
    """
    #define outputfiles
    pfamfile   = "{}.uniprot".format(PFAMID)
    intermediatefile = "{}.protnucl".format(PFAMID)
    full_length_nucleotides = "{}.fna".format(PFAMID)
    trimmed_nucleotides = "{}_trimmed.fna".format(PFAMID)

    # Download the PFAM Alignemnt file in Stockholm Format
    try:
        pfamsite   = "http://pfam.xfam.org/family/{}/alignment/uniprot".format(PFAMID)
        r = requests.get(pfamsite, stream=True)
        with open(pfamfile, 'wb') as fd:
            for chunk in r.iter_content(1000):
                fd.write(chunk)
    except:
        raise Error("Could not Download the PFAM File")

    # parse the file to obtain UniProt Accession IDs, start, and stop locations
    accessions = []
    for idx, line in enumerate(open(pfamfile,'r')):
        if line[0] not in ["#", "/"]:
            name, *stuff = line.split()
            try:
                accession, protrange = name.split("/")
                accession = accession.split(".")[0]
                start, stop = protrange.split("-")
                accessions.append((accession, start, stop))
            except:
                raise ValueError("Problem Parsing PFAM file at{}".format(name))

    # trimmed accessions
    ENA_protein_accessions = [acc for (acc, start, stop) in accessions]

    # write out intermediate file
    with open(intermediatefile,'w') as f:
        for (accession, start, stop) in accessions:
            f.write("{}\t{}\t{}\n".format(accession, start,stop))

    # obtain mappings of ENA files and mrege with protein data to get a four column
    # dataframe
    try:
        uniprot_ENA_mappings = uni.map(ENA_protein_accessions, f='ACC', t='EMBL')
        ENAaccessions = {k:v.pop() for k,v in uniprot_ENA_mappings.items()}
        df1 = pd.DataFrame(accessions, columns=["protein","start","stop"])
        df1 = df1.set_index("protein")
        df2 = pd.Series(ENAaccessions, name="Nucleotide")
        df3 = df1.join(df2)
        df3 = df3.dropna().reset_index()
        df3 = df3.rename(columns={"index":"Protein"} )
        df3.to_csv(intermediatefile)
    except:
        raise ValueError("Uniprot Mapping Failure")

    # Obtain Accessions from ENA
    try:
        ENAaccessions = [val.split(".")[0] for val in ENAaccessions.values()]
        getAccessionsfromENA(accessions=ENAaccessions, outfile=full_length_nucleotides)
    except:
        raise ValueError("Error in Obtaining Accessions from NCBI")

    # Trim the Nucleotides
    trimSequences(fna=full_length_nucleotides , translationfile=intermediatefile,
                  outfile=trimmed_nucleotides, PFAMID=PFAMID)

    #delete PFAM file if neccessary
    if not keepfiles:
        os.remove(pfamfile)
        os.remove(full_length_nucleotides)
        os.remove(intermediatefile)


def trimSequences(fna, translationfile, outfile, PFAMID):
    """
    Take an FNA file from ENA and a textfile with
    protein, start, stop, anc accesison numbers
    and trim the fata entries accordingly.
    """

    df = pd.read_csv(translationfile)
    print(df)

    nuclist = defaultdict(list)
    for idx, row in df.iterrows():
        nucl  = row["Nucleotide"]
        start = row["start"] * 3
        stop  = (row["stop"] * 3) + 2
        nuclist[nucl].append({"start":start, "stop":stop})

    #print(nuclist)
    recs = SeqIO.parse(open(fna,'r'),'fasta')

    outrecs = []
    for rec in recs:
        ena, accession, acc_ver = rec.id.split("|")
        try:
            locs = nuclist[acc_ver]
            for idx,loc in enumerate(locs):
                start = locs[0]["start"]
                stop  = locs[0]["stop"]
                rec.id = "{}|{}|{}".format(rec.id, PFAMID, idx)
                rec.seq = rec.seq[start:stop]
                outrecs.append(rec)
        except:
            pass

    print(outrecs[:5])
    SeqIO.write(outrecs, outfile, "fasta")



def getAccessionsfromENA(accessions, outfile, chunksize = 250):
    """
    Take a list of accessions that may be in the thousands,
    break it into smaller chunks and pull down the data from the
     European Nucleotide database using their REST api.
    """
    baseurl = "http://www.ebi.ac.uk/ena/data/view/"
    fasta =  "&display=fasta&download=fasta"

    #chunkedlist = [accession[i:i+chunksize] for i in range(0, len(accessions), chunksize)]
    chunkedaccessions = [tuple(accessions[x:y]) for (x, y) in [(x, x+chunksize) for x in range(0, len(accessions), chunksize)]]

    with open(outfile,'w') as f:
        for chunk in chunkedaccessions:
            names = ",".join(chunk)
            url = baseurl + names + fasta
            r = requests.get(url)
            if not r.status_code == 200:
                raise Exception('Bad response from Accessions {}'.format(names))

            else:
                print("Obtained Chunked Accessions; writing to file")
                f.write(r.text)

        print("All Accesions written!")



#getUniprotIDsfromPFAM("PF01867", outfile="PF01867.fna")
#trimSequences("PF01867.fna","PF01867.protnucl", outfile="testout.fna", PFAMID = "PF01867")