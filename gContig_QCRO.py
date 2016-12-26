import re
import linecache
import sys
import swalign

# This is for QCRO only

# Sample input:
#Best Structure= [35, 17, 28, 0, 10, 46, 31, 3, 11, 41, 47, 24, 38, 19, 1, 23, 6, 36, 27, 4, 32, 37, 2, 44, 12, 29, 34, 14, 43, 5, 45, 39, 9, 33, 18, 22, 40, 13, 25, 7, 15, 16, 20, 8, 30, 42, 21, 26]

# value for scoring
scoring = swalign.NucleotideScoringMatrix(1,-3)
sw = swalign.LocalAlignment(scoring)


# example of usage
# python gContig.py record.file frag.file

# input record file
infile = sys.argv[1]

# input DNA Fragments file
infile2 = sys.argv[2]

fo = open(infile, "r")

lines = fo.readlines()
#print lines[2]

# Separate lines[2]

p=[]
newstr= lines[0].replace(",", "")
newstr= newstr.replace("[", "")
newstr= newstr.replace("]", "")

x=newstr.split()
for s in x:
	#print s
	if s.isdigit():
		p.append(int(s))

# p is the list of name of DNA fragments that is arranged by the code
print p

fo.close()

#-----------
#Open fragment file

#fi = open("frag_x60189_4.dat", "r")

# write to file the newly arranged fragments
splitName = infile2.split("_", 1)

w2fArr = "Arr_{0}".format(splitName[1])
write2file = open(w2fArr, "w")

fragfile = infile2

for es in p:

	l=linecache.getline(fragfile, ((es*2)+1) )
	l = l.strip()
	print>>write2file, l

	l=linecache.getline(fragfile, ((es*2)+1+1) )
	l = l.strip()
	print>>write2file, l

write2file.close()

print p 

# open susun.file
#RArr = open(w2fArr, "r")

# generate info

# value for scoring
scoring = swalign.NucleotideScoringMatrix(1,-3)
sw = swalign.LocalAlignment(scoring, -2)


seq_db = {}
klist = []
fo = open(w2fArr, "r")
for l in fo:

	if re.match('^\>', l):
		name = l.split()
		key = name[0]
		klist.append(key)
		#print key

	elif re.match('^(a|c|g|t|A|C|G|T)', l):
		#print l
		seq_db[key]=l


	elif re.match(r'^\s*$', l):
		continue

print "klist : ",klist

infoFileWR = "info_{0}".format(splitName[1])
write2file2 = open(infoFileWR, "w")

for i in range(len(klist)-1):
	
	#print seq_db[p[i]], seq_db[p[i+1]], p[i], p[i+1]
	a = sw.align(seq_db[klist[i]], seq_db[klist[i+1]], klist[i], klist[i+1]).dump()
	print >> write2file2, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9]
	

a = sw.align(seq_db[klist[-1]], seq_db[klist[0]], klist[-1], klist[0]).dump()
#print type(a)
print >> write2file2, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9]

write2file2.close()
fo.close()
# Finish write info file
#------------------------------------------------------------------------------------

# Make Layout

fiInfo = open(infoFileWR, "r")
fragInfo = fiInfo.readlines()

main_contig = []
offset = 0

for i in range( len(fragInfo) ) :


	fragInfo[i].strip()
	
	if i<1:
		continue

	yyy = fragInfo[i].split()
	#print yyy
	rc= yyy[0]
	qname=yyy[1]
	rname=yyy[2]
	qstart=int(yyy[3])
	qend= int(yyy[4])
	len_q= int(yyy[5])
	rstart= int(yyy[6])
	rend= int(yyy[7])
	len_r= int(yyy[8])
	score = int(yyy[9])

	#print rc, qname, rname, qstart, qend, len_q, rstart, rend, len_r, score
	
	#set current fragment
	seq = seq_db[qname]

	if i==0:
		main_contig=seq_db[rname]

	print main_contig

	# to set the offset value to merge the next sequence
	offset=offset + rstart

	
	if len_r > len_q:
		# replace
		tempitem = main_contig[(offset+(qend-qstart)+1):]
		del main_contig[offset:]
		main_contig.extend(seq[qstart:])
		main_contig.extend(tempitem)

	else:
		del main_contig[offset:]
		main_contig.extend(seq[qstart:])

	

	if len(seq[:qstart]) > len(main_contig[:offset]):
		pos = len(seq[:qstart]) - len(main_contig[:offset])
		y = seq[:pos]
		main_contig[0:0] = y
		#offset = 0

	
	offset = offset - len(seq[:qstart])

	# if rstrand is reverse
	if rstart>rend:
		# rname is reverse compliment
		currentOrientation = reverse_Orientation[currentOrientation]

# out of for loop
j = "";
Wgenome = "".join(main_contig)
Wgenome = "".join(Wgenome.split())

# write main contig to file
newTitle = infile.split(".")
newTitle2 = newTitle[0].split("c")
mContigFile = "{0} Draft gene QCRO.fa".format(newTitle2[1])
fconWrite = open(mContigFile, "w")
print >>fconWrite, ">Draft gene {0} QCRO".format(newTitle2[1])
print >> fconWrite, Wgenome

fconWrite.close()
fiInfo.close()