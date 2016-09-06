# Orthogroup_QC.py
# input[1] file (ex: OrthologousGroup.txt) from OrthoFinder
# input[2] minimum number of genes per group
# output: OrthologousGroups_QC.txt file in same directory as original

#import libraries
import sys

#make function OG_QC
def OG_QC(filename, number):
	out = open(filename+".QC", 'w')
	Ortho = open(filename, "r")
	count=0
	for line in Ortho:
		if len(line.split()) > number:
			out.write(line)
			count=count+1
	print "kept {} groups".format(count)
	Ortho.close()
	out.close()


# to run on comand line as python Orthogroup_QC.py
#This tells python to actually run the following code instead of just importing libraries and defining functions
if __name__ == "__main__":
	#read in script inputs
	file = sys.argv[1]
	num = int(sys.argv[2])

	print "Output file is "+ file + ".QC"
	print "keep groups with {} or more members".format(num)

	#run the function
	OG_QC(file, num)
