import pandas as pd
import glob
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_files", help='input_files')
parser.add_argument("output_files",help='output_files')
args=parser.parse_args()
input = args.input_files
output = open(args.output_files,'w')

filez = glob.glob(input + "*.cnt")
print(filez[1])
t1 = pd.read_csv(filez[0], header=0, sep='\t')
tout = t1.iloc[:,0]
for f in filez:
	t1= pd.read_csv(f, header=0, sep='\t')
	tout= pd.concat([tout, t1.iloc[:,6]], axis=1)

tout.to_csv(output)
