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
import pysam

#This is the class library that is used with other small RNA analysis tools in this suite.

class SmallRNABam(object):
	def __init__(self, object):
		self.bam		= pysam.AlignmentFile(object, 'rb')
		self.clusters		= dict()
		self.filteredclusters	= dict()
		self.hairpins		= dict()
		self.matures		= dict()
		self.depths		= dict()
		self.hairpindepths	= dict()
		self.folds		= dict()
		self.hairpindepths	= dict()
		self.hairpinpositions	= dict()
		self.clusterdepths	= dict()
		self.clusterpositions	= dict()
		self.mdepths		= dict()
	def getClusters(self):
		name			= ''
		previousposition	= 0
		count			= 1
		start			= 0
		end			= 0
		columncount		= 0
		checkcount		= 0
		for column in self.bam.pileup():
			columncount	+= 1
		for column in self.bam.pileup():
			checkcount	+= 1
			if checkcount == columncount:
				self.clusters[count] = name, start, end
			if column.pos != previousposition + 1:
				if name != '':
					if previousposition != 0:
						self.clusters[count]	= name, start, end
						start	= column.pos
						count 	+= 1
			else:
				end		= column.pos
			previousposition	= column.pos
			name			= column.reference_name
	def getMatures(self, cutoff):
		mcount			= 0
		for key, value in self.clusters.iteritems():
			sequences		= dict()
			for read in self.bam.fetch(value[0], value[1], value[2]):
				if read.flag == 16 or read.flag == 272:
					strand	= '-'
				else:
					strand	= '+'
				positions	= [read.reference_name, read.reference_start, read.reference_end, strand, read.seq]
				positionsstring	= str(read.reference_start)+str(read.reference_end)+str(strand)
				if positionsstring not in sequences:
					sequences[positionsstring]	= [1]
					sequences[positionsstring].extend(positions)
				else:
					sequences[positionsstring][0]	+= 1
			sortedmatures	= sorted(sequences.items(), key=lambda x: x[1], reverse=True)
			for i in range(0, len(sortedmatures)):
				item	= sortedmatures[i]
				if int(item[1][0]) >= int(cutoff):
					mcount += 1
					if key not in self.filteredclusters:
						self.filteredclusters[key]	= value
						item	= list(item)
						item[1].append('m_'+str(mcount))
						self.matures[key]		= list(item[1:])
						self.mdepths['m_'+str(mcount)]	= item[1][0]
					else:
						item	= list(item)
						item[1].append('m_'+str(mcount))
						self.matures[key].extend(list(item[1:]))
						self.mdepths['m_'+str(mcount)]	= item[1][0]

	def extractHairpins(self, reference):
		fasta	= open(reference, 'r')
		fasta	= SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
		for key in sorted(self.filteredclusters.iterkeys()):
			start	= int(self.filteredclusters[key][1])
			end	= int(self.filteredclusters[key][2])
			if (end - start) < 100:
				end	= end + 60
				start	= start - 60
			if start <= 0:
				start	= 1
			if end >= len(fasta[self.filteredclusters[key][0]].seq):
				end	= len(fasta[self.filteredclusters[0]].seq)
			seq	= fasta[self.filteredclusters[key][0]].seq[start:end]
			if self.matures[key][0][4] == '-':
					seq	= seq.reverse_complement()
			self.hairpins[key]	= self.filteredclusters[key][0], start, end, self.matures[key][0][4], str(seq)
	def getHairpinPositions(self):
		fasta		= open('hairpins.fa', 'r')
		hairpins	= dict()
		for line in fasta.readlines():
			if re.match('>', line):
				header	= line[1:]
				header	= header.split('_')
				header	= header[0]
			else:
				if not re.match('\.|\(|\)', line):
					line	= re.sub('U', 'T', line)
					if int(header) not in hairpins:
						hairpins[int(header)]	= [line[:-1]]
					else:
						hairpins[int(header)].append(line[:-1])
		fasta.close()
		for key, value in self.hairpins.iteritems():
			extract	= value[-1]
			if key in hairpins:
				self.folds[key]	= []
				for item in hairpins[key]:
					start	= extract.find(str(item))
					endi	= start + len(str(item))
					if value[3] == '-':
						end	= value[2] - start
						start	= value[2] - endi
					else:
						start	= value[1] + start
						end	= value[1] + endi
					self.folds[key].extend([value[0], start, end, value[3], item])
	def getHairpinDepths(self, filteredfolds):
		for item in filteredfolds:
			for column in self.bam.pileup(item[1], item[2], item[3]):
				if item[0] in self.hairpindepths:
					self.hairpindepths[item[0]].append(column.nsegments)
					self.hairpinpositions[item[0]].append(column.pos)
				else:
					self.hairpindepths[item[0]]	= [column.nsegments]
					self.hairpinpositions[item[0]]	= [column.pos]
	def getClusterDepths(self):
		for key, value in self.clusters.iteritems():
			for column in self.bam.pileup(value[0], value[1], value[2]):
				if key in self.clusterdepths:
					self.clusterdepths[key].append(column.nsegments)
					self.clusterpositions[key].append(column.pos)
				else:
					self.clusterdepths[key]	= [column.nsegments]
					self.clusterpositions[key]	= [column.pos]
