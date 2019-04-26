from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Input_file", help='path to FASTA file with reads')

args=parser.parse_args()
Input_file = open(args.Input_file,'r')

Result = open('CpGisland.txt', 'w')
header = 'ID\tstart\tend\tstrand\tC%+G%\tObs/Exp\tlen\n'
Result.write(header)

def Cproc(seq):
	count = 0
	for i in seq:
		if i=='C':
			count = count + 1
	return float(count*100/len(seq))

def Gproc(seq):
	count = 0
	for i in seq:
		if i=='G':
			count = count + 1
	return float(count*100/len(seq))

def CGcount(seq):
	count=0
	for i in range(0,len(seq)-1):
		if seq[i]=='C' and seq[i+1]=='G':
			count = count+1
	return float(count)

def CGproc(seq):
	Cpro = Cproc(seq)
	Gpro = Gproc(seq)
	return float(Cpro + Gpro)

def Ccount(seq):
	count=0	
	for i in seq:
		if i == 'C':
			count=count+1
	return float(count)

def Gcount(seq):
	count=0	
	for i in seq:
		if i == 'G':
			count=count+1
	return float(count)

def CGtopattern(seq):
	C = Ccount(seq)
	G = Gcount(seq)
	expect = float(C*G/len(seq))
	CG = float(CGcount(seq))
	if expect != 0:
		return round(float(CG/expect),2)

def check(seq):
	if CGproc(seq)>50 and CGtopattern(seq)>0.6:
		return True


def CpGwindow(seq,window,ID,start,end,strand):
	start = 0
	end = 0
	wynik = []
	if(strand == '-'):
		seq = seq[::-1] # change direction to 5' -> 3'
	for i in range(0,len(seq)-window):
		if start==0: start = i
		if(check(seq[i:i+window])):
			end = i+window
		else:
			if(end!=0):	
				wynik.append([start,end])
			start = 0
			end = 0

	for i in wynik:
		if(i[1]-i[0]>200):
			Result.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(ID, i[0],i[1],strand,CGproc(seq[i[0]:i[1]]),CGtopattern(seq[i[0]:i[1]]),i[1]-i[0],'\n'))
			
	

for seq_record in SeqIO.parse(Input_file, "fasta"):

	seq = seq_record.seq
	ID = seq_record.id
	desc=seq_record.description
	start = desc.find('start:')
	end = desc.find('end')
	strand = desc.find('strand')
	start= desc[start+6:end-1]
	end= desc[end+4:strand-1]
	strand = desc[strand+7:strand+8]
	CpGwindow(seq,100,ID,start,end,strand)
	
	
