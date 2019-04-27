from Bio import SeqIO
import statistics
from statistics import median
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Input_file", help='path to FASTA file with reads')

args=parser.parse_args()
Input_file = open(args.Input_file,'r')

Result = open('CpGislands.txt', 'w')
header =  'ID\tstart\tend\tstrand\tC%+G%\tObs/Exp\tlen\n'
Result.write(header)


def CpGScanner(seq):
	cpgposition = []
	for i in range(0,len(seq)-1):
		if seq[i] == 'C' and seq[i+1] == 'G':
			cpgposition.append(i)
	return cpgposition


def distance(seq):
	cpgposition = CpGScanner(seq)
	dis = []
	for i in range(0,len(cpgposition)-1):
		dis.append(cpgposition[i+1] - cpgposition[i]-1)
	return int(median(dis))


def IslandMaker(start,end,island):
	if len(island) == 0:
		island.append([start,end])
	elif start-island[-1][1] <= 100:
		island[-1][1] = end
	else:
		island.append([start,end])
	return island


def CpGisland(seq):
	start = 0
	end = 0
	island = []
	dis = distance(seq)
	CpG = CpGScanner(seq)
	for i in range(0,len(CpG)-1):
		if CpG[i+1] - CpG[i] <= dis + 1:
			if start == 0: start = CpG[i]
			end = CpG[i+1] + 2
		elif end!=0:
			IslandMaker(start,end,island)
			start = 0
			end = 0
		if i == len(CpG) - 2:
			IslandMaker(start,end,island)
	return island	


def Cproc(seq):
	count = 0
	for i in seq:
		if i == 'C':
			count = count + 1
	return float(count*100/len(seq))


def Gproc(seq):
	count = 0
	for i in seq:
		if i == 'G':
			count = count + 1
	return float(count*100/len(seq))


def CGcount(seq):
	count = 0
	for i in range(0,len(seq)-1):
		if seq[i] == 'C' and seq[i+1] == 'G':
			count = count + 1
	return float(count)


def CGproc(seq):
	Cpro = Cproc(seq)
	Gpro = Gproc(seq)
	return float(Cpro + Gpro)


def Ccount(seq):
	count = 0	
	for i in seq:
		if i == 'C':
			count = count + 1
	return float(count)


def Gcount(seq):
	count = 0	
	for i in seq:
		if i == 'G':
			count = count + 1
	return float(count)


def CGtopattern(seq):
	C = Ccount(seq)
	G = Gcount(seq)
	expect = float(C*G/len(seq))
	CG = float(CGcount(seq))
	if expect != 0:
		return round(float(CG/expect),2)


def check(seq):
	if CGproc(seq) >= 50 and CGtopattern(seq) >= 0.6:
		return True		


def CpGCluster(seq,ID,start,end,strand):
	if(strand == '-'):
		seq = seq[::-1]
	island = CpGisland(seq)
	for i in island:
		if i[1] - i[0] > 200:
			if check(seq[i[0]:i[1]]):
				Result.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(ID, i[0]+1+start,i[1]+1+end,strand,CGproc(seq[i[0]:i[1]]),CGtopattern(seq[i[0]:i[1]]),i[1]-i[0],'\n'))


for seq_record in SeqIO.parse(Input_file, "fasta"):
	seq = seq_record.seq
	desc=seq_record.description
	desc = desc.strip()
	splitted_desc = desc.split(';')
	ID = splitted_desc[0].split(':')
	start = splitted_desc[1].split(':')
	start = start[1].split(':')	
	end = splitted_desc[2].split(':')
	end = end[1].split(':')
	strand = splitted_desc[3].split(':')
	strand = strand[1].split(':')
	strand = ''.join(strand)
	ID = ''.join(ID)
	start = ''.join(start)
	end = ''.join(end)

	CpGCluster(seq,ID,int(start),int(end),strand)
	
	
