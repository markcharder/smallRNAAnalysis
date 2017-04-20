#!/usr/bin/env python

import re
import argparse
import smallRNALib
import numpy
import glob
import os
import sys


#Simple Python script for gleaning information on small RNA sequencing data from fastq files.

parser	= argparse.ArgumentParser()
parser.add_argument(	'-a',
			'--fastq',
			dest='fastq',
			metavar='FASTQ',
			help='Input fastq file or directory containing multiple files.'	)
options	= parser.parse_args()

if not options.fastq:
	parser.error('Please provide fastq input file.')

if os.path.isdir(options.fastq):
	files	= os.listdir(options.fastq)
	print 'Sample', 'Size', 'A', 'G', 'C', 'T', 'total'
	if options.counts:
		fastq	= smallRNALib.SmallRNAFastq(options.fastq + '/' + file)
		fastq.getBiases()
		for key in fastq.starts.iterkeys():
			total = numpy.sum(fastq.starts[key]) 
			print file, key, fastq.starts[key][0], fastq.starts[key][1], fastq.starts[key][2], fastq.starts[key][3], total
else:
	fastq	= smallRNALib.SmallRNAFastq(options.fastq)
	print 'Sample', 'Size', 'A', 'G', 'C', 'T', 'total'
	sys.stdout.flush()
	fastq.getBiases()
	for key in fastq.starts.iterkeys():
		total = numpy.sum(fastq.starts[key]) 
		print str(os.path.basename(options.fastq)), key, fastq.starts[key][0], fastq.starts[key][1], fastq.starts[key][2], fastq.starts[key][3], total
		sys.stdout.flush()
