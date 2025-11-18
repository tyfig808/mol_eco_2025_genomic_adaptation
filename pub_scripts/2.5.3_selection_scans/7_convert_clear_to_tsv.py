# script to convert all clear output to tsv file
import sys
import pandas as pd
import os
import glob
import re

# create directory
directory = '~/mol_eco_2024/2.5.3_sel_scans'

# create df list to save to
#DF_list= list()

# loop over those with this pattern
for cf in glob.glob(directory + 'clear*.out'):
    
	# split based on _  and take the last line with the treat 
	result = re.split("_", cf)
	r=result[2]

	# split on the . and take treat
	t = re.split("[.]", r)
	treat = t[0]

	# read in clear output
	input=pd.read_pickle(cf)

	# add treatment column
	#input = input.insert(-1,'treatment', next_tave_val)

	# write out to tsv
	out_d = directory + "clear_" + treat + ".tsv"
	print(out_d)
	input.to_csv(out_d, sep="\t")

	# add data frames to list
	#DF_list.append(input)

#clear_lhb.tsv

# make total and write out
#result = pd.concat(frames)