## test.fasta
simple fastafile used to test location

## fq1.fq and fq2.fq
There are simple fastq files that can be used for testing concatenation
    seqtk sample -s100 AD13AK03_S1_L001_R1_001.fastq.gz 40 > fq1.fq
    seqtk sample -s100 AD13AK03_S1_L001_R2_001.fastq.gz 40 > fq2.fq
    

## trimmed_small.fq barcodes_small.fq and MappingFile.txt
These are three files for testing the parallel split function.
Its 1000 fastq sequences froma recent test run. The trimmed file is the fastq file while
th barcode file has the barcodes, the information for which is in MappingFile.txt
 
parallel_split_libraries_fastq should be able to split these files along their barcodes
using multiple cpus.