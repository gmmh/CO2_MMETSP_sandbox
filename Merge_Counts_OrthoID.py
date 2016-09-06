#Merge_Counts_OrthoID.py
# This script will take geneID.txt files and combine with _cds_counts.txt files

import re
import os
import sys

def merge_counts_orthoID(geneID_path, countfile_path):
	#Open count file for this taxa
	counts = open(countfile_path, 'r')
	cl = counts.readlines()
	header = cl.pop(0) #store the first line as a header
	
	#open gene ID file
	gID = open(geneID_path, 'r')
	
	#create outfile for writing
	outfile_path = "".join(re.findall("(\S*).txt",countfile_path)) + "_orthoID.txt"
	print "writing file {}".format(outfile_path)
	outfile = open(outfile_path, 'w')

	#make new header for outfile
	finalhead=['PeptideID','cds','Orthogroup']+header.split()
	outfile.write("\t".join(finalhead))
	outfile.write("\n")

	#loop through gene ID file, check that cds ID is same in counts file, append the count info onto the gene ID info and write to outfile
	for line in gID:
		y = line.split()
		cds= y[1]
		for n in range(len(cl)):
			if cds == ''.join(re.findall('[|]([0-9]+_[0-9]+)',cl[n])):
				x = cl.pop(n).split()
				#x = cl[n].split()
				break
		out = y + x
		outfile.write("\t".join(out))
		outfile.write("\n")
	counts.close()
	gID.close()
	outfile.close()

# This is the business end of the script where the functions get called.
if __name__ == "__main__":
	#Get inputs (1) gene ID directory, (2) count directory
	geneID_directory = sys.argv[1]
	count_directory = sys.argv[2]

	#pull out filenames from geneID_directory
	geneID_files = re.findall("\S*_geneIDs.txt", " ".join(os.listdir(geneID_directory)))
	print "Found geneID files: {}".format(geneID_files)

	#pull out filenames from count_directory
	count_files = re.findall("\S*_cds_counts.txt", " ".join(os.listdir(count_directory)))
	print "Found count files: {}".format(count_files)

	#Check if length of file lists are same, run function on each pair.
	if len(geneID_files) == len(count_files): #if there are same number of each
		for i in range(len(geneID_files)):
			gID_path = "".join(geneID_directory + geneID_files[i])
			c_path = "".join(count_directory + count_files[i])
			merge_counts_orthoID(gID_path, c_path)
	else:
		print "Error, geneID files not matched by count files"




