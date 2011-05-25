import os, sys, re

# read in arguments from command line
if len(sys.argv) != 3:
    print "need directory names"
    print "Usage: dats2DATs.py <source_dir> <target_dir>"
source_dir_name = sys.argv[1]
target_dir_name = sys.argv[2]

# this function opens a .dat formatted file and returns a string with .DAT file format
def dat2DAT(file_name, num_cols, num_rows):
	num_data_rows_seen = 0
	DAT_content = ""
	file = open(file_name,'r')

	for line in file.readlines() :
		if re.match('[\d\-\.\s]+\n', line) and not re.match('[\s\-]*\n', line):
			# If the line contains only digits, negative signs, dots, and white space
			# and does not contain ony dashes and whitespace convert the line to .DAT format
			DAT_line = re.sub('\s+\.(\s)',r'\tinf\1', line)
			DAT_line = re.sub('\s+\-','\t-', DAT_line)
			DAT_line = re.sub('\s+(\d)',r'\t\1', DAT_line)
			DAT_line = re.sub('\s+(\.\d)',r'\t\1', DAT_line)
			DAT_line = DAT_line.lstrip().rstrip()+"\n"
			DAT_content += DAT_line

	return DAT_content

for file_prefix in ["dangle", "int11", "int21", "int22", "loop", "miscloop", "sint2", "sint4", "stack", "tstackh", "tstacki", "tstackm"]:
	file_name = source_dir_name+"/"+file_prefix+".dat"
	output_file_path = target_dir_name+"/"+file_prefix+".DAT"
	output_file = open(output_file_path, 'w')
	output_file.write(dat2DAT(file_name, 16, 16))
	#print dat2DAT(file_name, 16, 16)

file_name = source_dir_name+"/tloop.dat"
output_file_path = target_dir_name+"/tloop.DAT"
DAT_content = ""
file = open(file_name,'r')

for line in file.readlines() :
	if re.match('[\d\-\.\sACGU]+\n', line) and not re.match('[\s\-]*\n', line):
		# If the line contains only digits, negative signs, dots, and white space or ACGU
		# and does not contain ony dashes and whitespace convert the line to .DAT format
		DAT_line = re.sub('\s+\.(\s)',r'\tinf\1', line)
		DAT_line = re.sub('\s+\-','\t-', DAT_line)
		DAT_line = re.sub('\s+(\d)',r'\t\1', DAT_line)
		DAT_line = re.sub('\s+(\.\d)',r'\t\1', DAT_line)
		DAT_line = DAT_line.lstrip().rstrip()+"\n"
		DAT_content += DAT_line
output_file = open(output_file_path, 'w')
output_file.write(DAT_content)
	
