#! /usr/bin/env python
# python selecthumanmouseonly.py
# author Vanessa Vermeirssen
#This script processes the bam files of reads mapping to the human genome and the mouse genome in parallel.
#It creates new "human preferred over mouse" and "mouse preferred over human" read list files.
#files as input: samplename_sroutH/Aligned.sortedByCoord.out.bam and samplename_sroutM/Aligned.sortedByCoord.out.bam
#output files: samplename_sroutH/HnotM_strict_querynames.txt (reads ONLY mapping to human - not further processed by exRNAxeno),
#samplename_sroutH/HnotM_querynames.txt, samplename_sroutM/MnotH_strict_querynames.txt (reads ONLY mapping to mouse - not further processed by exRNAxeno),
#samplename_sroutM/MnotH_querynames.txt

import pysam
from collections import defaultdict
import sys
import re
import os

##generate pairs of properly paired reads - assumes only uniquely mapped reads!
def read_pair_generator(bam):
    # """
    # Generate read pairs in a BAM file.
    # Reads are added to read_dict until a pair is found.
    # """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam:
        keep = True
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        #remove multiple aligned reads
        for entry in read.tags:
            if ('NH' in entry and entry[1] > 1):
                keep = False
        if keep:
            qname = read.query_name
            if qname not in read_dict:
                if read.is_read1:
                    read_dict[qname][0] = read
                else:
                    read_dict[qname][1] = read
            else:
                if read.is_read1:
                    yield read, read_dict[qname][1]
                else:
                    yield read_dict[qname][0], read
                del read_dict[qname]

#input files
fileH=sys.argv[1]
fileM=sys.argv[2]
#output files
file_out1=sys.argv[3]
file_out2=sys.argv[4]
file_out3=sys.argv[5]
file_out4=sys.argv[6]
samfileH = pysam.AlignmentFile(fileH, "rb")
samfileM = pysam.AlignmentFile(fileM, "rb")
sample=re.search("(RNA\d{6}_srout\w)",fileH).group()
[sample1, sample2] = sample.split("_srout")

##calculate the edit distance for both BAM/SAM files
#The edit distance of a paired read is defined as the sum of NM tags (Number of differences,
#mismatches plus inserted and deleted bases, between the sequence and reference) and
#soft-clipping from both read pairs.
#p = 0 #number of properly paired mapped reads from primary alignment
distH={}
distHspikes={} #put Sequin and ERCC spikes in separate dictionary
for read1, read2 in read_pair_generator(samfileH): #AlignedSegment object
    #p += 1
    readName1 = read1.query_name
    readName2 = read2.query_name
    readRef1 = read1.reference_name
    readRef2 = read2.reference_name
    editDist = 0
    tag1 = 0
    tag2 = 0
    for item in read1.cigartuples: #cigartuples 4 - soft-clipping; NM-tag 10 if present from aligner
        if (item[0] == 4):
            editDist += item[1]
    try:
        tag1 = read1.get_tag("NM")
    except KeyError:
        tag1 = 0
    editDist = editDist + tag1
    for item in read2.cigartuples: #cigartuples 4 - soft-clipping
        if (item[0] == 4):
            editDist += item[1]
    try:
        tag2 = read2.get_tag("NM")
    except KeyError:
        tag2 = 0
    editDist = editDist + tag2
    #Edit distance is defined as sum of NM tags and soft-clipping from both read pairs
    #print(readName1, read1.get_tag("NM"), read2.get_tag("NM"), editDist)
    if (re.search("ERCC", readRef1) or re.search("chrIS", readRef1)):
        distHspikes[readName1] = 1
    else:
        distH[readName1]= editDist
    if (re.search("ERCC", readRef2) or re.search("chrIS", readRef2)):
        distHspikes[readName2] = 1
    else:
        distH[readName2]= editDist
    #print(distH[readName1], distH[readName2])

#pm = 0 #number of properly paired mapped reads from primary alignment
distM={}
distMspikes={}
for read1, read2 in read_pair_generator(samfileM):
    #pm += 1
    readName1 = read1.query_name
    readName2 = read2.query_name
    readRef1 = read1.reference_name
    readRef2 = read2.reference_name
    editDist = 0
    tag1 = 0
    tag2 = 0
    for item in read1.cigartuples:
        if (item[0] == 4):
            editDist += item[1]
    try:
        tag1 = read1.get_tag("NM")
    except KeyError:
        tag1 = 0
    editDist = editDist + tag1
    for item in read2.cigartuples:
        if (item[0] == 4):
            editDist += item[1]
    try:
        tag2 = read2.get_tag("NM")
    except KeyError:
        tag2 = 0
    editDist = editDist + tag2
    if (re.search("ERCC", readRef1) or re.search("chrIS", readRef1)):
        distMspikes[readName1] = 1
    else:
        distM[readName1]= editDist
    if (re.search("ERCC", readRef2) or re.search("chrIS", readRef2)):
        distMspikes[readName2] = 1
    else:
        distM[readName2]= editDist

#print("There are {} number of properly paired mapped reads to human".format(p))
#print("There are {} number of properly paired mapped reads to mouse".format(pm))

samfileH.close()
samfileM.close()

nameHn=[] #query_name not mouse according to some edit distance
nameMn=[] #query_name not human according to some edit distance

m=0
n=0
k=0
l=0
##compare human and mouse edit distances - human
for i in (set(distH.keys()) - set(distM.keys())):
    nameHs.append(i)

for i in (set(distH.keys()) & set(distM.keys())): #so if a read pair maps equally well to the human and mouse ref genome, it is discarded
    if (distH[i] < distM[i]):
        nameHn.append(i)

#add spike reads back
for i in distHspikes.keys():
    nameHs.append(i)

nameHnM = nameHn + nameHs
m=len(nameHs)
n=len(nameHnM)

##compare human and mouse edit distances - mouse
for i in (set(distM.keys()) - set(distH.keys())):
    nameMs.append(i)

for i in (set(distH.keys()) & set(distM.keys())): #so if a read maps equally well to the human and mouse ref genome, it is discarded
    if (distM[i] < distH[i]):
        nameMn.append(i)

#add spike reads back
for i in distMspikes.keys():
    nameMs.append(i)

nameMnH = nameMn + nameMs
k=len(nameMs)
l=len(nameMnH)
#print(len(nameHs))
#print(nameHs)
#print(len(nameHnM))
#print(sample1, p, pm, m, n, k, l)

##writing output to files
with open(file_out1, 'w') as f:
    for item in nameHs:
        f.write(item+"\n")

with open(file_out2, 'w') as f:
    for item in nameHnM:
        f.write(item+"\n")

with open(file_out3, 'w') as f:
    for item in nameMs:
        f.write(item+"\n")

with open(file_out4, 'w') as f:
    for item in nameMnH:
        f.write(item+"\n")
