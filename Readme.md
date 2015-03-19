### Tools for Working with Qiime.


###install
```[python]
pip install --editable .
```


#### fastqconcat
**fastqconcat** will take two fastq files and concatenate each record set keeping only a certian amount of sequence and
optionally reverse/complementing the reverse sequence/quality score

```[bash]
fastqconcat --forward_fastq fastq_f.fastq.gz \
            --reverse_fastq fastq_r.fastq.gz \
            --outfile out.fq \
            --keep_left 250 \
            --keep_right 175 \
            --ncpus 2
```

#### parallel_concat
**parallel_concat** will take two fastq files and concatenate each record. The method is different as it splits
 the files using system `split` and `cat` functions rather than using BioPython's SeqIO iterators.
```[bash]
parallel_concat --forward_fastq fastq_f.fastq.gz \
            --reverse_fastq fastq_r.fastq.gz \
            --outfile out.fq \
            --keep_left 250 \
            --keep_right 175 \
            --ncpus 2
            --splitsize 1000000 #howmany lines to split on
```

#### parallel_split_library_fastq
**parallel_split_library_fastq** wraps Qiime's split_libraries_fastq.py program. This script will split the files,
process them in parallel, and aggregate the results. Because it can run in parallel it us faster - MUCH faster.
However this currently only supports one fastq file and one barcode file.

```[bash]
parallel_concat --fastq reads.fastq \
            --barcode_fastq barcodes.fastq \
            --outfile seqs.fna \
            --mappingfile Mapping.txt \ #Qiime Mapping File
            --barcodetype 8 \
            --qual_cutoff 19 \
            --keep_left 250 \
            --splitsize 100000 \ # number of lines to split on (must be a multiple of 4)
            --keep_right 175 \
            --ncpus 30 #if you got em....
```
