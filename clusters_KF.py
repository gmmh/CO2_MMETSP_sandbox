import sys  


#this is a script that harriet and i wrote to turn the output from MCL clustering into
#a transcript-to-gene-map that can be input into RSEM

#first you make a dictionary, call it OG
def getOG(filename):
	OG = {}
	prefix = "OG_"
	#filename = "/Users/kyle_frischkorn/Desktop/test.txt"
	handle = open(filename,"r")
	for i,line in enumerate(handle):
		if i%1000==0:
			print "I'm working bitch, I've done ",i
		trans_list = line.split("\t")
		name = prefix+str(i)
		for item in trans_list:
			trans_id = item#.split("#")[0][:-1]
			if name in OG.keys():
				OG[name].append(trans_id)
			else:
				OG[name]=[trans_id]
	return OG
	
				
def makefile(OG,filename):
	outfile = open(filename+".out","w")
	for key in OG.keys():
		for item in OG[key]:
			outfile.write("\t".join([key,item]))
			outfile.write("\n")
	outfile.close()

#the steps below makes it so that you can run the script as a program using the command:
# python cluster.py INPUT	
if __name__ == "__main__":
	fname=sys.argv[1]
	
	print "I'm taking MCL output file and turning it into Transcript-to-gene-map for RSEM mapping :) "
	outog=getOG(fname)
	makefile(outog,fname)
	