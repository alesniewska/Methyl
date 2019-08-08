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
t1 = pd.read_table(filez[1], header=1)
tout = t1.iloc[:,0]
for f in filez:
	t1= pd.read_table(f, header=1)
	tout= pd.concat([tout, t1.iloc[:,6]], axis=1)
tout.to_csv(output)
