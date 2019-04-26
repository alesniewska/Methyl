import argparse

parser = argparse.ArgumentParser()
parser.add_argument("CT_conversion", help='CT_conversion file - from bismark')
parser.add_argument("GA_conversion", help='GA_conversion file - from bismark')
parser.add_argument("Count_Reads_file", help='File with counted mapped reads')

args=parser.parse_args()

CT_conversion  = open(args.CT_conversion,'r')
GA_conversion  = open(args.GA_conversion,'r')
Openfile  = open(args.Count_Reads_file,'r').readlines()

firstLine=Openfile.pop(0)

final = []
CT = []
GA = []
seq = []


for line in Openfile[1:]:

	line = line.strip()
	splitted = line.split('\t')
	geneID = splitted[0].split(';')
	start = splitted[2].split(';')
	end = splitted[3].split(';')
	strand = splitted[4].split(';')
	for i in range(0,len(start)):
		for a,b,c,d in zip(geneID,start[i],end[i],strand[i]):
			final.append((a,start[i],end[i],d))
final.sort(key=lambda x:x[1], reverse=True)

for line in CT_conversion:
	line = line.strip()
	line = line.strip()
	if(line.startswith('>')):
		continue
	else:
		for char in line:
			CT.append(char)
for line in GA_conversion:
	line = line.strip()
	if(line.startswith('>')):
		continue
	else:
		for char in line:
			GA.append(char)


for elem in final:
	info = '>{};start:{};end:{};strand;{}{}'.format(elem[0],elem[1],elem[2],elem[3],'\n')
	
	if elem[3] == '-':
		sequence = GA[int(elem[1]):int(elem[2])]
		sequence = ''.join(sequence).strip()

	else:
		sequence = CT[int(elem[1]):int(elem[2])]
		sequence = ''.join(sequence).strip()
		
	seq.append((info,sequence.strip()))


open('reads.fasta', 'w').write('\n'.join('%s%s' % x for x in seq))

