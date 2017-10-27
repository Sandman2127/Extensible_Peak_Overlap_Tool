#!/usr/bin/env python
from __future__ import division
import math
import sys

"""
	usage: python Gene_map_overlap.py <TAIR10_genes_only.bed> <TF_narrowpeakfile> <promoter or gene> <bp_upstream>
								by: Dean Sanders
								  Oct 27,2017
"""


# arg =
# chr1	3631	5899		0	+	gene	AT1G01010

#arg2=
#chr3	910184	910385	3:910284	1000	.	2416.0	0.00	999.00	100

arg = sys.argv[1]
arg2 = sys.argv[2]
arg3 = str(sys.argv[3]).strip()  # possible options map_all , promoter_only , gene_only
try:
	arg4 = int(sys.argv[4])
except IndexError:
	pass

print "argument total is:",len(sys.argv)


def readict(arg):
	with open(arg, 'rb') as fin:
		rows = (line.split('\t') for line in fin)
		dictA = {row[0]+row[1]+row[6]:row[1:] for row in rows}
		return dictA

def all_overlap(dictA,dictB,bp):
	counter = 0
	for keys in dictB:
		counter += 1
		chrom_B = str(keys[:4]).strip()
		peakstart = int(dictB[keys][0])
		peakend = int(dictB[keys][1])
		for key in dictA :
			if dictA[key][4].strip() == "+" :
				gene_name = dictA[key][6].strip()
				strand = dictA[key][4].strip()
				chrom_A = str(key[:4]).strip()
				upstream = (int(dictA[key][0]) - bp)
				gene_start = int(dictA[key][0])
				gene_end = int(dictA[key][1])
				downstream = (gene_end + bp)
				if (chrom_A == chrom_B) & (peakstart <= upstream) & (peakend >= upstream): # overlap left side region
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"overlap_promoter_front",gene_name
				elif (chrom_A == chrom_B) & (peakstart >= upstream ) & (peakend <= gene_start): # within bp region
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"within_promoter",gene_name
				elif (chrom_A == chrom_B) & (peakstart <= gene_start) & (peakend >= gene_start): # within bp region into gene
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"overlap_gene_border_TSS",gene_name
				elif (chrom_A == chrom_B) & (peakstart >= gene_start) & (peakend <= gene_end): # within the gene
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"within_gene",gene_name
				elif (chrom_A == chrom_B) & (peakstart >= gene_end) & (peakend <= downstream): # downstream of the gene by arg4 bp
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"downstream",gene_name
				else:
					pass

			if dictA[key][4].strip() == "-" :
				gene_name = dictA[key][6].strip()
				strand = dictA[key][4].strip()
				chrom_A = str(key[:4]).strip()
				upstream = (int(dictA[key][1]) + bp)
				gene_start = int(dictA[key][1])  # reversed due to change in direction
				gene_end = int(dictA[key][0])   # same
				downstream = (gene_end - bp)
				if (chrom_A == chrom_B) & (peakstart <= upstream) & (peakend >= upstream): # overlap right side region
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"overlap_promoter_front",gene_name
				elif (chrom_A == chrom_B) & (peakstart >= gene_start) & (peakend <= upstream): # within bp region
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"within_promoter",gene_name
				elif (chrom_A == chrom_B) & (peakstart <= gene_start) & (peakend >= gene_start): # within bp region running into gene
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"overlap_gene_border_TSS",gene_name
				elif (chrom_A == chrom_B) & (peakstart >= gene_end) & (peakend <= gene_start): # within the gene
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"within_gene",gene_name
				elif (chrom_A == chrom_B) & (peakstart >= downstream) & (peakend <= gene_end): # downstream of the gene by arg4 bp
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"downstream",gene_name
				else:
					pass

def promoter_only(dictA,dictB,bp):
	counter = 0
	for keys in dictB:
		counter += 1
		chrom_B = str(keys[:4]).strip()
		peakstart = int(dictB[keys][0])
		peakend = int(dictB[keys][1])
		for key in dictA :
			if dictA[key][4].strip() == "+" :
				gene_name = dictA[key][6].strip()
				strand = dictA[key][4].strip()
				chrom_A = str(key[:4]).strip()
				upstream = (int(dictA[key][0]) - bp)
				gene_start = int(dictA[key][0])
				gene_end = int(dictA[key][1])
				if (chrom_A == chrom_B) & (peakstart <= upstream) & (peakend >= upstream): # overlap left side region
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"overlap_promoter_front",gene_name
				elif (chrom_A == chrom_B) & (peakstart >= upstream ) & (peakend <= gene_start): # within bp region
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"within_promoter",gene_name
				elif (chrom_A == chrom_B) & (peakstart <= gene_start) & (peakend >= gene_start): # within bp region into gene
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"overlap_gene_border_TSS",gene_name
#
			if dictA[key][4].strip() == "-" :
				gene_name = dictA[key][6].strip()
				strand = dictA[key][4].strip()
				chrom_A = str(key[:4]).strip()
				upstream = (int(dictA[key][1]) + bp)
				gene_start = int(dictA[key][1])  # reversed due to change in direction
				gene_end = int(dictA[key][0])   # same
				if (chrom_A == chrom_B) & (peakstart <= upstream) & (peakend >= upstream): # overlap right side region
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"overlap_promoter_front",gene_name
				elif (chrom_A == chrom_B) & (peakstart >= gene_start) & (peakend <= upstream): # within bp region
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"within_promoter",gene_name
				elif (chrom_A == chrom_B) & (peakstart <= gene_start) & (peakend >= gene_start): # within bp region running into gene
					print chrom_A,chrom_B,upstream,gene_start,gene_end,peakstart,peakend,strand,"overlap_gene_border_TSS",gene_name

def genebody_map(dictA,dictB):
	counter = 0
	for keys in dictB:
		counter += 1
		chrom_B = str(keys[:4]).strip()
		peakstart = int(dictB[keys][0])
		peakend = int(dictB[keys][1])
		for key in dictA :
			if dictA[key][4].strip() == "+" :
				gene_name = dictA[key][6].strip()
				strand = dictA[key][4].strip()
				chrom_A = str(key[:4]).strip()
				gene_start = int(dictA[key][0])
				gene_end = int(dictA[key][1])
				if (chrom_A == chrom_B) & (peakstart >= gene_start) & (peakend <= gene_end): # within the gene
					print chrom_A,chrom_B,gene_start,gene_end,peakstart,peakend,strand,"within_gene",gene_name
#
			if dictA[key][4].strip() == "-" :
				gene_name = dictA[key][6].strip()
				strand = dictA[key][4].strip()
				chrom_A = str(key[:4]).strip()
				gene_start = int(dictA[key][1])  # reversed due to change in direction
				gene_end = int(dictA[key][0])   # same
				if (chrom_A == chrom_B) & (peakstart >= gene_end) & (peakend <= gene_start): # within the gene
					print chrom_A,chrom_B,gene_start,gene_end,peakstart,peakend,strand,"within_gene",gene_name

if __name__ == "__main__":
	if len(sys.argv) == 5 :
		if arg3 == "map_all" :
			all_overlap(readict(arg),readict(arg2),arg4)
		elif arg3 == "promoter_only":
			promoter_only(readict(arg),readict(arg2),arg4)
		if arg3 == "genebody_only":
			print "bp argument is irrelevant with genebody_only option, mapping to genes only"
			genebody_map(readict(arg),readict(arg2))
		else:
			pass
	if len(sys.argv) == 4 :
		if arg3 == "genebody_only":
			genebody_map(readict(arg),readict(arg2))
		else:
			pass
	#except IndexError:
		#print "usage: python Gene_map_overlap.py <TAIR10_genes_only.bed> <TF_narrowpeakfile> <overlap_type> <bp_upstream> \n overlap_types = map_all , promoter_only, genebody_only"
