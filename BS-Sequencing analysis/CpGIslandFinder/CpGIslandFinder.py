from Bio import SeqIO
import statistics
from statistics import median
import argparse


parser = argparse.ArgumentParser(
  description='''CpGIslandFinder is a program to find CpG island HOW TO USE: python3 CpGIslandFinder.py <path to FASTA file> <optional arguments>  ''')


parser.add_argument("fasta_file", help='path to FASTA file')
parser.add_argument('-min_island_length', type=int, default=500, help='Minimal island length, default is 500 (int type)')
parser.add_argument('-obs_to_exp', type=float, default=0.6, help='Obs/Exp value, default is 0.6 (float type)')
parser.add_argument('-Cproc_plus_Gproc', type=float, default=50.0, help='C procent + G procent in sequence, default is 50.0 (float type)')


args=parser.parse_args()
min_island_length = args.min_island_length
obs_to_exp = args.obs_to_exp
Cproc_plus_Gproc = args.Cproc_plus_Gproc

fasta_file = open(args.fasta_file,'r')
Result = open('CpGislands.gtf', 'w')
#header =  'chromosome\tstart\tend\tstrand\tC%+G%\tObs/Exp\tCG_count\tlen\n'
#Result.write(header)


#distance beetwen CpG in sequence
def Distance(seq):
	cpgposition = CpGScanner(seq)
	dis = []
	for i in range(0,len(cpgposition)-1):
		dis.append(cpgposition[i+1] - cpgposition[i] -1)
	if len(dis)!=0:
		return int(median(dis))
	else: return 0

# join CpGislands if distance beetwen islands next to each other is lower than 100
def IslandMaker(start,end,island):
	if len(island) == 0:
		island.append([start,end])
	elif start-island[-1][1] < 100:
		island[-1][1] = end
	else:
		island.append([start,end])
	return island

#finding CpGislands if distance beetwen CpG is lower than distance calculated by Distance(seq) 
def CpGisland(seq):
	start = 0
	end = 0
	island = []
	dis = Distance(seq)
	CpG = CpGScanner(seq)
	for i in range(0,len(CpG)-1):
		if CpG[i+1] - CpG[i] -1 < dis:
			if start == 0: start = CpG[i]
			end = CpG[i+1] + 2
		elif end!=0:
			IslandMaker(start,end,island)
			start = 0
			end = 0
		if i == len(CpG) - 2 and end!=0:
			IslandMaker(start,end,island)
	return island	

#finding CpG positions
def CpGScanner(seq):
	cpgposition = []
	for i in range(0,len(seq)-1):
		if seq[i] == 'C' and seq[i+1] == 'G':
			cpgposition.append(i)
	return cpgposition

def Cproc(seq):
	seq = seq.upper()
	count = 0
	for i in seq:
		if i == 'C':
			count = count + 1
	return float(count*100/len(seq))


def Gproc(seq):
	seq = seq.upper()
	count = 0
	for i in seq:
		if i == 'G':
			count = count + 1
	return float(count*100/len(seq))

#calculate the number of CpG in sequence
def CGcount(seq):
	count = 0
	for i in range(0,len(seq)-1):
		if seq[i] == 'C' and seq[i+1] == 'G':
			count = count + 1
	return float(count)


def CGproc(seq):
	Cpro = Cproc(seq)
	Gpro = Gproc(seq)
	return round(float(Cpro + Gpro),2)


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

# CpG count / expected , expected = (C count * G count / length of sequence)
def CGtopattern(seq):
	C = Ccount(seq)
	G = Gcount(seq)
	expect = float(C*G/len(seq))
	CG = float(CGcount(seq))
	if expect != 0:
		return round(float(CG/expect),2)

# checking is C% + G+ equal or bigger than Cproc_plus_Gproc (default is 50.0) 
# and CpG count / expected is equal or bigger than obs_to_exp value (default is 0.6)
def check(seq):
	if CGproc(seq) >= Cproc_plus_Gproc and CGtopattern(seq) >= obs_to_exp:
		return True		

def CpGCluster(seq,desc,strand):
	island = CpGisland(seq) #finding islands
	for i in island:
		if i[1] - i[0] > min_island_length: # checking by min_length
			if check(seq[i[0]:i[1]]): # is check output is true
					Result.write('CpGisland_{}_{}-{}\t{}\tCpGisland\t{}\t{}\t.\t{}\t.\tLength {}; C%+G% {}; Obs/Exp {};\n'.format(desc,i[0]+1,i[1]+1,desc,i[0]+1,i[1]+1,strand,i[1]-i[0],CGproc(seq[i[0]:i[1]]),CGtopattern(seq[i[0]:i[1]]))) #island is writing in outputfile!
				
#island is writing in outputfile!!

#chromosome = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
for seq_record in SeqIO.parse(fasta_file, "fasta"):
	seq = seq_record.seq
	seq = seq.upper()
	desc=seq_record.description
	CpGCluster(seq,desc,'+')
	seq = seq.reverse_complement()
	CpGCluster(seq,desc,'-')

	
