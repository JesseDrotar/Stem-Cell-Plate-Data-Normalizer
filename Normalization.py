#!/usr/bin/python
#This program reads in .csv data from 384 well plates read for Luciferase signal by an automated plate reader
#The program reads in both a experimental plate (BGplate) and a control plate (CTGplate)
#The values are normalized respectively and then a ratio is calculated that represents the sample's overall signal fold induction

import sys
import csv
import collections
from collections import OrderedDict

BGplate = sys.argv[1]
CTGplate = sys.argv[2]


Bright_Glow_List = []
CTG_List = []


with open(BGplate, 'r') as BG:

	reader = csv.reader(BG, delimiter=',')

	for i in reader:
		Bright_Glow_List.append(i)

with open(CTGplate, 'r') as CTG:
	
	reader = csv.reader(CTG, delimiter=',')
	
	for i in reader:
		CTG_List.append(i)

#screens for non int values in the file	
BGtotal = 0
BGcounter = 0
converted_Bright_Glow_List = []

for i in Bright_Glow_List:
	for j in i:
		j = int(j)
		converted_Bright_Glow_List.append(j)
		BGtotal = BGtotal + j
		BGcounter = BGcounter + 1

#Calculates the average value for the bright glow plate
BGaverage = BGtotal/BGcounter

#Creates a filter value
BGplate_filter = BGaverage * 10

#Screens data for any outliers above the filter limit
filtered_Bright_Glow_List = []
for i in converted_Bright_Glow_List:
	if (i > BGplate_filter):
		pass
	else:
		filtered_Bright_Glow_List.append(i)

#The last column of the plate is used as a value for normalization
#This totals the values of the last column
BGbackground_luminescent = filtered_Bright_Glow_List[24] + filtered_Bright_Glow_List[48] + filtered_Bright_Glow_List[72] + filtered_Bright_Glow_List[96] + filtered_Bright_Glow_List[120] + filtered_Bright_Glow_List[144] + filtered_Bright_Glow_List[168] + filtered_Bright_Glow_List[192] + filtered_Bright_Glow_List[216] + filtered_Bright_Glow_List[240] + filtered_Bright_Glow_List[264] + filtered_Bright_Glow_List[288] + filtered_Bright_Glow_List[312] + filtered_Bright_Glow_List[336]

BGbackground_average = BGbackground_luminescent/ 14

normalized_Bright_Glow_List = []

#Normalizes the values
for i in filtered_Bright_Glow_List:
	normalize_val = i - BGbackground_average
	normalized_Bright_Glow_List.append(normalize_val)

#Repeats the same process for the Control plates
CTGtotal = 0
CTGcounter = 0
converted_CTG_List = []
for i in CTG_List:
	for j in i:
		j = int(j)
		converted_CTG_List.append(j)
		CTGtotal = CTGtotal + j
		CTGcounter = CTGcounter + 1
	
CTGaverage = CTGtotal/CTGcounter

CTGplate_filter = CTGaverage * 10

filtered_CTG_List = []
for i in converted_CTG_List:
	if (i > CTGplate_filter):
		pass
	else:
		filtered_CTG_List.append(i)
		
CTGbackground_luminescent = filtered_CTG_List[24] + filtered_CTG_List[48] + filtered_CTG_List[72] + filtered_CTG_List[96] + filtered_CTG_List[120] + filtered_CTG_List[144] + filtered_CTG_List[168] + filtered_CTG_List[192] + filtered_CTG_List[216] + filtered_CTG_List[240] + filtered_CTG_List[264] + filtered_CTG_List[288] + filtered_CTG_List[312] + filtered_CTG_List[336]

CTGbackground_average = CTGbackground_luminescent/14

normalized_CTG_List = []

for i in filtered_CTG_List:
	normalize_val = i - CTGbackground_average
	normalized_CTG_List.append(normalize_val)
	

BG_CTG_Ratio_List = []

#Calculates the ratio of BG to CTG
for i, j in zip(normalized_Bright_Glow_List, normalized_CTG_List):
	BG_CTG = (i / j)
	BG_CTG_Ratio_List.append(BG_CTG)

#Totals the DMSO values for the first column of the plate
DMSO_Total = BG_CTG_Ratio_List[25] + BG_CTG_Ratio_List[49] + BG_CTG_Ratio_List[73] + BG_CTG_Ratio_List[97] + BG_CTG_Ratio_List[121] + BG_CTG_Ratio_List[145] + BG_CTG_Ratio_List[169] + BG_CTG_Ratio_List[193] + BG_CTG_Ratio_List[217] + BG_CTG_Ratio_List[241] + BG_CTG_Ratio_List[265] + BG_CTG_Ratio_List[289] + BG_CTG_Ratio_List[313] + BG_CTG_Ratio_List[337]

DMSO_Average = DMSO_Total / 14

Signal_Fold_Induction_List = []

for i in BG_CTG_Ratio_List:
	i = i/DMSO_Average
	Signal_Fold_Induction_List.append(i)
	
Signal_Fold_Induction_List = [Signal_Fold_Induction_List[x:x+24] for x in range(0, len(Signal_Fold_Induction_List), 24)]

#Creates a new user designated file that will contain the new normalized data in the same 384 well format.
file_name = input("Enter the final file name as a .csv file: ")

with open(file_name, 'w') as output:
	writer = csv.writer(output)
	
	for i in range(0,16) :
		writer.writerow(Signal_Fold_Induction_List[i])
