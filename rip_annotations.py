#Annotate_orthoID_geneID_index.py

#import libraries
import re
from ftplib import FTP #import ftp tools
import os #need for listing file names from orthoID_geneID_index and annot directory
import sys #need for taking inputs such as pep and annot directories

def rip_annotations(ID, save_dir):
	#access ftp and rip pfam annotations
	ftp= FTP('ftp.imicrobe.us') #set ftp server
	ftp.login() #log in
	ftp.cwd('camera/combined_assemblies') #change to main working directory
	ripdir= ID+"/annot"
	ftp.cwd(ripdir)
	#write to README file in working directory
	savefile= save_dir +ID+"_pfam.txt"
	print "Saving file {}".format(savefile)
	ftp.retrbinary('RETR pfam.gff3', open(savefile, 'wb').write)
	print "Done"
	ftp.quit()

#This annotation function does not yet work, just use the rip annotations function for now
def add_annotations(ID, index_dir, annot_dir):
	geneID_file = index_dir + ID + "_geneIDs.txt"
	geneIDs = open(geneID_file, "r")
	annot_file = annot_dir + ID + "_pfam.txt"
	annotations = open(annot_file, "r")
	pfams = annotations.readlines()
	outfile_path = annot_dir + ID + "_geneIDs_annot.txt"
	outfile = open(outfile_path, "w")
	for line in geneIDs:
		y = line.split()
		cds = y[1]
		x= ["NA"]
		for n in range(len(pfams)):
			if re.search(cds, pfams[n]):
				if cds == ''.join(re.findall('[|]([0-9]+_[0-9]+)',pfams[n])):
					x = pfams.pop(n).split()
		out = y + x
		outfile.write("\t".join(out))
		outfile.write("\n")
	geneIDs.close()
	pfams.close()
	outfile.close()

#run the script:
if __name__ == "__main__":
	index_dir = sys.argv[1]
	annot_dir = sys.argv[2]
	
	#make a list of all geneID index files by taxa name
	allfiles = os.listdir(index_dir)
	taxa = re.findall( "(\S*)_geneIDs.txt"," ".join(allfiles))
	print taxa
	
	#make list of all taxa with annotation files with pfam in the name
	annotfiles = os.listdir(annot_dir)
	taxa_w_pfams = re.findall("(\S*)_pfam", " ".join(annotfiles))
	print taxa_w_pfams
	
	for taxon in taxa:
		if re.search(taxon, " ".join(taxa_w_pfams)) == None:
			print "No pfam file for {}".format(taxon)
			try: 
				rip_annotations(taxon, annot_dir)
				print "Find annotations for {} in {}".format(taxon, annot_dir)
			except:
				print "Error for {}".format(taxon), sys.exc_info()[0]
				#break
		#add_annotations(taxon, index_dir, annot_dir)
