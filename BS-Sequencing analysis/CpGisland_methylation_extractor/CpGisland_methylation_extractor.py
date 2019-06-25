import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''CpG island methylation extractor''')
parser.add_argument("islands", help='path to islands file')
parser.add_argument("BedGraph", help='path to BedGraph file')
parser.add_argument("genom", help='path to genom file')

args=parser.parse_args()
islands = open(args.islands,'r')
BedGraph = open(args.BedGraph,'r')
genom = open(args.genom,'r')

Result = open("CpGisland_methylated.txt",'w')
Result.write("chromosome\tfeature\tstart\tend\tstrand\tC%+G%\tObs/Exp\tCGcount\tmethyl%\tgene\n")


#if methylation position point to CpG
def Methyl_position(seq,desc):
	methylation = []
	for line in BedGraph:
		splitted = line.split('\t')
		if splitted[0] == desc:
			if seq[int(splitted[1])] == 'C' and seq[int(splitted[2])] == 'G':
				methylation.append([int(splitted[1]),splitted[3]]) #remember position and percentage of methylation
	return methylation

#number of CpG in island
def counter(seq):
	count = 0
	for i in range (0,len(seq)-1):
		if seq[i] == 'C' and seq[i+1] == 'G':
			count = count + 1
	return count

#count percentage of methylation for CpG islands
def island_methylation(seq,desc,strand):
	methylation = []
	methylation = Methyl_position(seq,desc)
	for island in islands:
		island = island.split('\t')
		if island[4]==strand:
			count=0		
      for i in range(0,len(methylation)):
				if int(methylation[i][0]) >= int(island[2]) and int(methylation[i][0]) <= int(island[3]):
				  count = count + float(methylation[i][1])/100
			methyl = round(float((count)/float(island[7])),5) 
			Result.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(island[0],island[1],island[2],island[3],island[4],island[5],island[6],island[7],methyl,island[8]))


for seq_record in SeqIO.parse(genom, "fasta"):
	seq = seq_record.seq
	seq = seq.upper()
	desc=seq_record.description
	island_methylation(seq,desc,'+')
	seq = seq.reverse_complement()
	island_methylation(seq,desc,'-')
