#!/usr/bin/env python
import argparse
import re

def get_args():
    parser= argparse.ArgumentParser()
   # parser = argparse.ArgumentParser(description="This code is used to identify duplicates in SAM file reads. On single end reads, and is assuming the UMIs have been added to the read.THis code takes into consideration all possile cigar strings (inclduing soft clipping), strandedness, single-end reads, and known UMI")
    parser.add_argument("-f", "--filename", help="Input filename", required=True)
    parser.add_argument("-u", "--umilist", help="Input umi list", required=True)
    parser.add_argument("-o", "--outfile", help="Output filename", required=False)

    return parser.parse_args()

 
args = get_args()
f=args.filename
u=args.umilist
o=args.outfile



###########
#functions#
###########

def is_minus_strand(bit_flag):
    '''this function will look at the bitwise flag and return direction of strand, if it is a minus strand it will return true, if it is plus strand it will return false '''
    return (bit_flag & 0x10) != 0

def adjust_position(line, cigar):
    '''This function calculates and returns the adjusted position, taking into account the original position, the number of matches (M), soft clipping (S), deletions (D), and skips (N) in the CIGAR string. '''
    position = int(line[3])  # Extract and convert the position field to an integer
    num_M = 0
    num_D = 0
    num_N = 0
    right_soft_clip = 0
    left_soft_clip=0
    bit_flag=int(line[1])
    #determine which strand its on
    strand=is_minus_strand(bit_flag)
    #regular expression to find the right handed S value
    pattern = r'(\d+)S$'
    left_pattern = r'(\d+)S(?=\d+M)'


    #search for the pattern in the cigar string 
    match = re.search(pattern, cigar)
    left_match=re.findall(left_pattern, cigar)
    # Split the CIGAR string into individual operations, creates a list of tuples, where each tuple contains the length and the operation type for each operation in the CIGAR string.
    cigar_parts = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    #length holds the length of the current operation, and operation holds the type of the current operation.
    # operation it is (e.g., 'M' for match, 'D' for deletion, 'N' for skip, 'S' for soft clip, etc.).
    #length is the number associated with cigar, operation is the letter assosciated with the cigar
    for length, operation in cigar_parts:
        length = int(length)
        if operation == 'M':
            num_M += length
        elif operation == 'D':
            num_D += length
        elif operation == 'N':
            num_N += length
        #minus strand so only count right hand soft clipping
        elif operation == 'S' and strand and match:
            #need this to add only S' from the right side of M
            right_soft_clip = int(match.group(1))
           
        #left strand so only count left soft clipping        
        elif operation == 'S' and not strand and left_match:
            #need this to add only S's from the left siide of M
            left_soft_clip = int(left_match[-1])
    #adjust position
    if strand:
        adjusted_position = position + num_M + int(right_soft_clip) + num_D + num_N
    else:
        adjusted_position= position + num_D + num_M + num_N - left_soft_clip

    return adjusted_position

def has_soft_clipping(line):
    '''This function checks if there is any soft clipping in the CIGAR string of a SAM line.'''
    cigar = line.split()[5]
    # Check if there is 'S' in the CIGAR string
    if 'S' in cigar:
        return True 
    return False

####################################################################################################################
#file will be sorted on command line by chromosome, samtools sort C1_SE_uniqAlign.sam -o sorted_C1_SE_uniqAlign.sam#
####################################################################################################################

###set to hold as a tuple non-duplicate values (chrome, strand, position, umi) non-duplicates will also be written to an outfile. Set will be emptied when the chromosomes no longer are equal, to help with memory usage###
alignment=set()


###bring umi's in and store in a list, for comparison, if there are unknown UMI's it will be written to unknown_umi_file###
umi_list=[]
with open(u,"r") as file:
   for line in file:
      line = line.strip("\n")
      umi_list.append(line)
#close(u)

#initialize the varaibles  want to keep chrome, strand, position, umi
chrome=None
cigar=None
qname=None
strand=None
position=None
umi=None
prev_chrome=None

#duplicate count 
duplicate=0 
#open files
with open(f,"r") as f, open(o,"w") as outfile, open("unknown_umi.txt", "w") as unknown_umi_file:
    for line in f:
        if line.startswith('@'): 
            outfile.write(line)
            continue
#alignment line, put the first SAM file read into the alignment set, this will be our comparison point for duplicates  
        line=line.split("\t")
#needs to hold the very first line in the outfile 
        bit_flag=int(line[1])
        #if its unmapped continue
        if (bit_flag & 4) == 4:
                continue
        strand=is_minus_strand(bit_flag)
        qname=line[0]
        cigar=line[5]
        chrome=line[2]
        umi=line[0][-8:]
        #save unknown umi's to a file 
        if umi not in umi_list:
            unknown_umi_file.write(str(line))
            continue
        if is_minus_strand:
            position = adjust_position(line, cigar)
        else:
            if has_soft_clipping:
              position = adjust_position(line, cigar)
            else:
                position = int(line[3])
        elements=(chrome, strand, position, umi)
        
# Check if the elements tuple is already in the alignment list. If present in the alignment set, then it is a duplicate, and we want to not include in our outfile so we will continue. Otherwise the elements will be added to the alignment set  
        if elements in alignment:
            duplicate +=1 
            continue
        else:
            alignment.add(elements)
            outfile.write("\t".join(line))

        if chrome != prev_chrome: 
            alignment.clear()
            prev_chrome = chrome
print(duplicate)
#Close all files

# close(f)
# close(unknown_umi_file)
# close(o)
        
            

## to check after, on command line 
# grep "^@" output.txt (64)
# grep -v "^@" output.txt | wc -l ()
            

