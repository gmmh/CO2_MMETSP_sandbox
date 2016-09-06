#!/usr/bin/python
#MMETSP_sample_import.py
#inputs: (1) taxonomic classification of interest, (2) mmetsp taxa file with path
#outputs: rips and writes peptide fasta files and count data for each taxa

#import libraries
import sys
from ftplib import FTP #import the ftp library
import re #import regular expression tools

#Get arguments from command line
t=sys.argv[1] #full or partial taxa name
f=sys.argv[2] #mmetsp taxa file with path

#Load the mmetsp taxa file
mtf=open(f,'r')

#loop through taxonomy file looking for species name for taxonomic name supplied
g=[] #make an empty list to store genus names
for line in mtf:
	if re.search(t,line): #if taxa name in line
		g= g+line.split('\t')[7:8]#pull out the 8th field should be genus, keeping as list

g=set(g) #keep only unique genus names
print g
#close the taxonomy file
mtf.close()

#Access the ftp server for MMETSP combined assemblies
ftp= FTP('ftp.imicrobe.us') #set home ftp server
ftp.login() #log in
ftp.cwd('camera/combined_assemblies') #change working directory

#List files and find files matching genus
files=ftp.nlst() #make a list of all files and directories in wd
delimiter=' '
all=delimiter.join(files)
for genus in g:
	string= genus+"\S*.pep.fa.gz"
	taxafiles=re.findall(string, all)
	print "{} files matching genus=".format(len(taxafiles))+genus
	print taxafiles
	if len(taxafiles) > 0:
		for file in taxafiles:
			command = "RETR "+file
			outfile = "../work/"+file
			ftp.retrbinary(command, open(outfile, 'wb').write)

ftp.close()

