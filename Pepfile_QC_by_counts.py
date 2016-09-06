#Pepfile_QC_by_counts.py
#input [1]: pep.fa file from MMETSP
#input [2]: count file from MMETSP
#input [3]: minimum number of counts to consider good contig (suggest 2-10)
#output: qc.fasta file with peptide seqs removed that have no counts

import sys
import re
import os

def good_count_IDs(countfile, min_number):
	file = open(countfile, 'r')
	header = file.readline() #skips first line and stores header
	good_list=[]
	for line in file:
		#Then pop off the zeroith element of the list (the ID)
		count_list = line.split()
		ID = count_list.pop(0)
		#then sum the remaining elements (counts)
		counts = [int(x) for x in count_list]
		if sum(counts) >= min_number:
			#keep the popped off first element of the list
			good_list.append(ID)
	file.close()
	return good_list

def mk_qc_pepfile(good_list, pepfile):
	file =  open(pepfile, 'r') #open peptide file
	print file
	taxa = re.findall("(\S*).pep.fa", pepfile) #store taxa name
	output = open("".join(taxa+[".qc.pep.fa"]), 'w') #create outfile with same name .qc.pep.fa
	print output
	for line in file:
		ID = re.findall("NCGR_PEP_ID=(\S*)", line) #find the peptide ID
		#print ID
		if "".join(ID) in good_list: #if the pepID is in the good list
			#print "Found it"
			output.write(line)
			nextline = file.next()
			output.write(nextline)
	output.close()
	file.close()

#Execute the code below:
if __name__ == "__main__":
	#pep_file = sys.argv[1]
	pep_dir = sys.argv[1]
	#count_file = sys.argv[2]
	count_dir = sys.argv[2]
	min_number = int(sys.argv[3])
	
	pep_files = os.listdir(pep_dir)
	pep_files = re.findall("\S*(?!qc).pep.fa", " ".join(pep_files))
	print "Found these peptide files: \n {}".format(pep_files)
	
	count_files = os.listdir(count_dir)
	count_files = re.findall("\S*_counts.txt", " ".join(count_files))
	print "Found these count files: \n {}".format(count_files)
	
	if len(count_files) == len(pep_files):
		for i in xrange(0,len(count_files)):
			count_file = count_files[i]
			file_path = count_dir + "".join(count_file)
			print "Loading counts file {}".format(file_path)
			good_list = good_count_IDs(file_path, min_number)
			print "There were {} contigs with at least {} counts".format(len(good_list), min_number)
			taxa = re.findall("(\S*)_cds_counts.txt", count_file)
			pep_file = pep_files[i]
			if re.search("".join(taxa), "".join(pep_file)):
				file_path = pep_dir + "".join(pep_file)
				print "Loading peptide file {}".format(file_path)
				print "Making new QC peptide file"
				mk_qc_pepfile(good_list, file_path)
			else:
				print "Pepfile {} does not match countfile {}".format(pep_file, count_file)

