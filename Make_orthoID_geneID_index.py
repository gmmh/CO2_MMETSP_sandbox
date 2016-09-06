# Make_orthoID_geneID_index.py
#input [1] directory for '.pep.fa' files
#input [2] path and filename for OrthologousGroups.txt or OrthologousGroups_QC.txt made by OrthoFinder for these files
#outputs: a file for each pep.fa containing a list of PEP_ID, cds_ID and OrthoID

#load libraries
import re
import sys
import os

#Define a function to rip IDs from pep.fa and orthologous groups files
def mk_index(pepfile, orthofile):
	#open up files:
	pep = open(pepfile,'r')
	print pep
	base= re.findall('(\S*).pep.fa', pepfile)
	#make outfile
	outfile = open(''.join(base)+'_geneIDs.txt','w')
	print outfile
	#open up Ortholog file
	Ortho = open(orthofile, 'r')
	print Ortho
	#loop through pep.fa file to pull out IDs
	for line in pep:
		Ortho.seek(0) #start at beginning of Ortho file
		OrID=['NA'] #Null is no OrthoGroup
		if re.search('^>\S*', line):
			#pull out CAMPEP and cds IDS
			y=line.split() #split the line by whitespace
			#Need to add taxaID to 'out'
			PEPID=re.findall(">(CAMPEP_\S*)", y[0])
			cdsID=re.findall("([0-9]+_[0-9]+)", y[1])
			#Search Ortholog file for CAMPEP ID and pull out any orthoIDs
			for nl in Ortho:
				if re.search(''.join(PEPID), nl):
					#print "Found it!"
					OrID=re.findall("(OG[0-9]+):",nl)
					#print OrID
					break
			#concatenate these ortho IDs onto 'out', else concatenate 'NA'
			#Write to output file
			out= PEPID+cdsID+OrID
			outfile.write("\t".join(out))
			outfile.write("\n")
	pep.close()
	outfile.close()
	Ortho.close()
	print "Done making index for " + ''.join(base)

#Execute the code below:
if __name__ == "__main__":
	#read in variables from inputs
	directory = sys.argv[1]
	orthofile = sys.argv[2]
	
	files = os.listdir(directory)
	#print files
	pepfiles = re.findall("\S*.pep.fa", " ".join(files))
	#print pepfiles
	
	for pepfile in pepfiles:
		print "Running file {}".format(pepfile)
		filepath = directory + pepfile
		#print filepath
		#run the function
		mk_index(filepath, orthofile)
		print "Done"
