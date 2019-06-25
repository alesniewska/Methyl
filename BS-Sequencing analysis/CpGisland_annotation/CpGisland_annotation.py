from Bio import SeqIO
import statistics
from statistics import median
import argparse

parser = argparse.ArgumentParser(description='''CpG island annotation''')
parser.add_argument("CpGislands_file", help='path to CpGislands file')
parser.add_argument("exon_gtf_file", help='path to exon gtf file')

args=parser.parse_args()
CpGislands_file = open(args.CpGislands_file,'r')
gtf_file= open(args.exon_gtf_file,'r')

Result = open('CpGislands_annotation.txt', 'w')
header =  'chromosome\tstart\tend\tstrand\tC%+G%\tObs/Exp\tlen\tgene_id\n'
Result.write(header)

#get name of gene in gtf file
def getGeneName(string):
	length = len(r'gene_name "')
	start = string.find(r'gene_name "')
	stop = string.find(r'"',start + length)
	return string[start + length:stop]


#reading gtf file and writing in list
def gtf_annotation():
	annotation = []
	for line in gtf_file:
		line = line.strip()
		splitted = line.split('\t')
		gene_name = getGeneName(line)
		desc = 'chr' + splitted[0] 
		annotation.append([desc,splitted[2],splitted[3],splitted[4],splitted[6],gene_name])
	return annotation


#deleting duplicates in result
def deduplication(result):
	deduplicated = []
	for elem in result:
		if elem not in deduplicated:
			deduplicated.append(elem)
	return deduplicated


annotation = tuple(gtf_annotation())
result = []


for line in CpGislands_file:
	line = line.strip()
	splitted = line.split('\t')
	chromosome = splitted[0]
	if splitted[0] != 'chromosome':
		start = splitted[1]
		end = splitted[2]
		strand = splitted[3]
		for elem in annotation:
			if elem[0] == chromosome and elem[4] == strand and int(start)>=int(elem[2]) and int(end)<=int(elem[3]):
				result.append([chromosome,elem[1],start,end,strand,splitted[4],splitted[5],splitted[6],elem[5]])

deduplicated = deduplication(result)

for elem in deduplicated:
     Result.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(elem[0],elem[1],elem[2],elem[3],elem[4],elem[5],elem[6],elem[7],elem[8]))
