#!/usr/bin/python
#rip_counts_annot_MMETSP
#inputs: directory with .pep.fa files from MMETSP, outputs count and annot files for each taxa

#import libraries
import sys
import os
from ftplib import FTP #import ftp library
import re #import regular expression tools

# taxa directory
t=sys.argv[1]
files= os.listdir(t)

#Pull out taxa and strian ID from .pep.fa files
delimiter=' '
all=delimiter.join(files)
taxa=re.findall('(\S*).pep.fa',all)
print taxa

#access ftp
ftp= FTP('ftp.imicrobe.us') #set ftp server
ftp.login() #log in
ftp.cwd('camera/combined_assemblies') #change to main working directory


#This looop pulls out all count data for each taxa and saves in count directory
for ID in taxa:
    #change to taxa directory/readcounts
    ripdir= ID+"/readcounts"
    ftp.cwd(ripdir)
    #write to README file in working directory
    savefile= t +"/counts/"+ID+"_cds_counts.txt"
    ftp.retrbinary('RETR cds.dat', open(savefile, 'wb').write)
    ftp.cwd("~/camera/combined_assemblies")
#close ftp connection
ftp.quit()