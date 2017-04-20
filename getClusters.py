#!/usr/bin/env python

import re
import argparse
import smallRNALib
import numpy
import glob
import os
import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

#Simple Python script for gleaning information on small RNA sequencing data from fastq files.
#This script requires that RNAFold and hhmmir are installed and available in the system path. This includes the 'fold_regions.pl' script
#that is distributed with hhmmir.

parser  = argparse.ArgumentParser()
parser.add_argument(    '-a',
                        '--bam',
                        dest='bam',
                        metavar='BAM',
                        help='Input small RNA sequencing bam alignment file.' )
parser.add_argument(	'-b',
			'--reference',
			dest='reference',
			metavar='REFERENCE',
			help='Reference sequence from which to extract hairpins.'	)
parser.add_argument(	'-c',
			'--cutoff',
			dest='cutoff',
			metavar='CUTOFF',
			default=25,
			help='Cutoff depth to filter read mappings by, defaults to 25.'	)
parser.add_argument(	'-d',
			'--prefix',
			dest='prefix',
			metavar='PREFIX',
			default='out',
			help='Prefix of output files, defaults to "out".'	)

options = parser.parse_args()
foldpath='/usr/local/bin/RNAfold' #Set this to installation path for RNAfold binary from the Vienna suite.
if not options.bam:
	parser.error("Please specify input bam.")
if not options.reference:
	parser.error("Please specify input reference genome.")

print 'Reading in bam file...'
file	= smallRNALib.SmallRNABam(options.bam)
print 'Done'
print 'Finding read clusters from bam file...'
file.getClusters()
print 'Done'
print 'Finding all read stacks above cutoff value...'
file.getMatures(options.cutoff)
print 'Done'
print 'Extracting sequence regions surrounding read clusters...'
file.extractHairpins(options.reference)
print 'Writing to file...'
outfile	= open(options.prefix + '.extracts.fa', 'w')
for key, value in file.hairpins.iteritems():
	string	= '>' + str(key) + '\n' + str(value[-1]) + '\n'
	outfile.write(string)
print 'Done'

print 'Running RNAfold...'
#Run RNAfold on regions extracted from read clusters.
os.system('perl hhmmir_v1.2/fold_regions.pl ' + foldpath + ' ' + options.prefix + '.extracts.fa')
os.system('java -jar hhmmir_v1.2/ExtractHairpins.jar folded.fa hairpins.fa >> removed.txt')
print 'Done'

print 'Finding genomic positions of hairpins relative to clusters...'
file.getHairpinPositions()
print 'Done'

print 'Filtering out overlapping hairpins...'
foldslist	= list()
for key, value in file.folds.iteritems():
	for i in xrange(0, len(value), 5):
		foldslist.append([key] + value[i:i+4])
foldslist.sort(key=lambda x: (x[1], x[2]))
previous		= list()
filteredfoldslist	= list()
for item in foldslist:
	if previous != []:
		if item[1] == previous[1]:
			if item[2] >= previous[3]:
				filteredfoldslist.append(item)
		else:
			filteredfoldslist.append(item)
	previous	= item
print 'Done'

print 'Identifying sequences above threshold that overlap hairpins...'
count	= 1
for i in range(0, len(filteredfoldslist)):
	for key, value in file.matures.iteritems():
		for item in value:
			if int(item[2]) >= int(filteredfoldslist[i][2]) and int(item[3]) <= int(filteredfoldslist[i][3]) and item[1] == filteredfoldslist[i][1]: 
				filteredfoldslist[i].append(item[-1])
		if re.search('m_', str(filteredfoldslist[i][-1])):
			filteredfoldslist[i].append(key)
	filteredfoldslist[i][0]	= count
	count += 1
print 'Done'

print 'Generating depth files for cluster and hairpin regions...'
file.getClusterDepths()
file.getHairpinDepths(filteredfoldslist)
print 'Done'

print 'Printing results to file...'
outfile	= open(options.prefix + '.clusters.txt', 'w')
for key, value in file.clusters.iteritems():
	string	= map(str, value)
	string	= '\t'.join(string) + '\t'
	string	= string + str(key) + '\n'
	outfile.write(string)
outfile.close()
outfile	= open(options.prefix + '.hairpins.txt', 'w')
for item in filteredfoldslist:
	if re.search('m_', str(item[-2])):
		string	= map(str, item[1:])
		string	= '\t'.join(string) + '\t'
		string	= string + str(item[0]) + '\n'
		outfile.write(string)
outfile.close()
outfile		= open(options.prefix + '.matures.txt', 'w')
outfastafile	= open(options.prefix + '.matures.fasta', 'w')
for key, value in file.matures.iteritems():
	for item in value:
		string	= map(str, item[1:])
		string	= '\t'.join(string) + '\t'
		string	= string + str(item[0]) + '\n'
		fastastring	= '>' + str(item[-1]) + '\n' + str(item[5]) + '\n'
		outfile.write(string)
		outfastafile.write(fastastring)
outfile.close()
outfile	= open(options.prefix + '.clusterdepths.txt', 'w')
for key, value in file.clusterdepths.iteritems():
	for i in range(0, len(value)):
		string	= str(key) + '\t' + str(value[i]) + '\t' + str(file.clusterpositions[key][i]) + '\n'
		outfile.write(string)
outfile.close()
outfile	= open(options.prefix + '.hairpindepths.txt', 'w')
for key, value in file.hairpindepths.iteritems():
	for i in range(0, len(value)):
		string	= str(key) + '\t' + str(value[i]) + '\t' + str(file.hairpinpositions[key][i]) + '\n'
		outfile.write(string)
outfile.close()
outfile	= open(options.prefix + '.maturecounts.txt', 'w')
for key, value in file.mdepths.iteritems():
	string	= str(key) + '\t' + str(value) + '\n'
	outfile.write(string)
outfile.close()
print 'Done'
