import argparse

parser = argparse.ArgumentParser()

parser.add_argument("input_files", help='input_files')
parser.add_argument("id",help="id")
parser.add_argument("output_files",help='output_files')

args=parser.parse_args()
input = open(args.input_files,'r')
id = args.id
output = open(args.output_files,'w')

header = 'Geneid\tchr\tstart\tend\tstrand\tlength\t{}\n'.format(id)

print(header)
output.write(header)
for line in input:
	output.write(line)

