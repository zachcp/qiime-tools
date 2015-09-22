from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import zip

import multiprocessing
from functools import partial

import click
from Bio import SeqIO


@click.command()
@click.option('--fastq', type=click.File('r'), prompt=True,help="name of the fastq forward file")
@click.option('--outdir', prompt=True, help="name of the output directory")
def split_fastq_by_name(fastq,outdir):
	"""
	split a fastqfile based on the names in the header

	"""
	fastqs = SeqIO.parse(fastq, "fastq")

	for fastq in fastqs:
		#The record name is {{name}}_{{number}}
		sample = fastq.name.split("_")[0]
		outlocation = "{}/{}.fastq".format(outdir,sample)
		qualscores = map(str,fastq.letter_annotations['phred_quality'])

		with open(outlocation,'wa') as f:
			f.write("@{}\n{}\n+\n{}".format(fastq.description,
											fastq.seq,
											" ".join(qualscores))
	print("Splitting Complete")