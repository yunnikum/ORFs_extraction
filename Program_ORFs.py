############################################################################
#
# Function:     filesearchandread()
#
# Programmer:   Mary Mills
#
# Class:        BIFS 617 9040 Advanced Bioinformatics Summer 2020
#
# Assignment:   Group Project
#
# Purpose:      Collect user input and check if file exists in working 
#               directory. Isolate headers and sequence data from file.
#               
#               
#
# Date:         July 29, 2020
#               
#
# Version:      0.0
#                          
#
# Input:        DNA sequence in a text file
#
# Output:       Properly formatted DNA sequence with bases all uppercase
#               and extra characters removed
#               Header stored
#
# References:   
#               
#
############################################################################

import os
import re

#Check that the file is in the working directory and read the sequence
def filesearchandread():
    fileinput = input("What is the file name you want to search? ")
    print()
    #If the file exists in the directory, begin to process it
    if os.path.exists(fileinput) == True:
        print("File located", "\n")
        #Keep file open and begin to read the lines
        with open(fileinput, 'r') as seqfile:
            seqdata = seqfile.read()
            #REGEX to store accession number
            accnum = re.findall('>.+]', seqdata, re.MULTILINE)
            #REGEX to find sequence data, ignoring case type and whitespace
            seq = re.findall('^[a-zA-Z\s]{3,100000000}', seqdata, re.MULTILINE)
            #begin loop to process any number of sequences in the file
            for i in range(0,len(accnum)):
                #print accession numbers and sequences
                print(accnum[i], '\n', seq[i])    
    else:
        print('That file does not exist. Please try again')


########################################################################################
#
#Function: read_input_file_name() 
#
#Programmer: Yekaterina Unnikumaran
#
#Class: BIFS 617 9040 Advanced Bioinformatics Summer 2020
#
# Purpose: Determine if file(s) is/are present in working directory.
#
# Date: July 25, 2020
#
# Version: 1.0
#
# Input: Calling FASTA file(s)
#
# Output: Print FASTA sequences
#
# References: www.python.org,www.pythonforbiologists.com, www.w3schools.com/python
#
########################################################################################

import re
import string

#create the function and open the file(s)
def read_input_file_name():
    with open('FASTA1_GroupProject.txt', 'rb') as f:
        fasta1 = f.read()

    fasta1 = [x.split('\n', 1) for x in fasta1.split('>')]
    fasta1 = [(x[0], ''.join(x[1].split())) for x in fasta1 if len(x) == 2]
        
    sequence = fasta1[103:] #removes first 103 characters which are accession number and all else in first line
                            #that is not a part of the sequence
    
    seq = str(sequence.replace('\n','')) #removes the rest of the spaces in between the lines so only characters
                                         #left are the sequence characters, converted to strings
    print(seq) #print sequence ready to be used
    print()
    
########################################################################################
# Function: read_input_ORF_length
#
# Programmer: Yekaterina Unnikumaran
#
# Class: BIFS 617 9040 Advanced Bioinformatics Summer 2020
#
# Assignment: Group Project
#
# Purpose: Prompt user if ORF < 50 characters in length
#
# Date: July 25, 2020
#
# Version: 1.0
#
# Input: Call for sequence to separate sequence into 50 character chunks
#
# Ouput:Prints sequence counted off in 50 character groups
#
# References: www.python.org, www.pythonforbiologists.com, www.stackoverflow.com, www.ib.bioninja.com
# 
#Background on ORFs: To identify an open reading frame:
#
#Locate a sequence corresponding to a start codon in order to determine the reading frame Ã¢â‚¬â€œ this will be ATG (sense strand)
#Read this sequence in base triplets until a stop codon is reached (TGA, TAG or TAA)
#The longer the sequence, the more significant the likelihood that the sequence corresponds to an open reading frame
#
########################################################################################
#import re
import re

#create function
def read_input_ORF_length():
    file = open('FASTA1_GroupProject.txt', 'r') #opens the file
    fasta1 = file.read() #reads the file

#Cleaning the sequence so it can be used
    sequence = fasta1[103:] #removes first 103 characters which are accession number and all else in first line
                            #that is not a part of the sequence
    
    seq = str(sequence.replace('\n','')) #removes the rest of the spaces in between the lines so only characters
                                         #left are the sequence characters, converted to strings
    print('Sequence 1:',seq) #print sequence ready to be used
    print()
    print('Number of bases in sequence:', len(seq)) #count bases in sequence before we chunk it into groups of 50
                                                    #for ORF reading purposes
    print()

#To scan through the sequence for ORFs create an index
    print('Bases separated in groups of 50 to find ORFs with the next function \n', [seq[i:i+50] for i in range(0, len(seq), 50)])
    return seq

    
    
############################################################################
#
# Function:     determine_reading_frames_and_find_ORFs()
#
# Programmer:   Ben Sparklin
#
# Class:        BIFS 617 9040 Advanced Bioinformatics Summer 2020
#
# Assignment:   Group Project
#
# Purpose:      Determine 6 reading frames and store each reading frame in a list.
#               Search each reading frame for any ORFs and save any ORFs to orf_list.
#               
#
# Date:         July 22, 2020 / BS
#               
#
# Version:      1.0
#               
#               
#
# Input:        Forward strand sequence
#
# Output:       orf_list
#
# References:   www.stackoverflow.com
#               www.python.org
#
############################################################################

def determine_reading_frames_and_find_ORFs():

    import re

    #save forward strand
    seq = read_input_ORF_length()
    minimum_orf_length = seq

    #determine reverse strand sequence
    #determine complement of forward sequence
    forward_strand_complement = seq.translate(str.maketrans({"A":"T","C":"G","T":"A","G":"C"}))
    #reverse forward strand complement
    reverse_strand = forward_strand_complement[::-1]

    #Create empty reading frame list. The list will hold the 6 reading frames as strings
    #for each sequence in the following order: +1,+2,+3,-1,-2,-3
    reading_frame_list = []

    #add +1 reading frame to reading frame list with complete forward strand
    reading_frame_list.append(seq)

    #add +2 reading frame to reading frame list by removing first base pair and last two base pairs of forward strand
    reading_frame_list.append(seq[1:len(seq)-2])

    #add +3 reading frame to reading frame list by removing first two pairs and last base pair of forward strand
    reading_frame_list.append(seq[2:len(seq)-1])

    #add -1 reading frame to reading frame list with complete reverse strand
    reading_frame_list.append(reverse_strand)

    #add -2 reading frame to reading frame list by removing first base pair and last two base pairs of reverse strand
    reading_frame_list.append(reverse_strand[1:len(reverse_strand)-2])

    #add -3 reading frame to reading frame list by removing first two pairs and last base pair of reverse strand
    reading_frame_list.append(reverse_strand[2:len(reverse_strand)-1])

    #Set search pattern to find each ORF.
    #Match must start with ATG, has a three chracter sequence consisting of A,T,G, or C for 0 or more times,
    #and end with TAA, TAG, or TGA
    pattern = r"(?=ATG)(?:[ATGC]{3})*?(?:TAA|TAG|TGA)"

    #In each of the 6 reading frames, search for an ORF. If an ORF is found,and meets the length requirement save it to orf_list
    orf_list = []
    for i in range(6):
         match = re.findall(pattern, reading_frame_list[i]) #search for all ORF patterns in reading frame
         #if match: #if  ORF(s) are found in a given reading frame, save the the result to orf_list
         #    for a in range(len(match)):
          #       if len(match[a]) >= minimum_orf_length:#if an ORF is greater than or equal the minimum ORF length, add it to orf_list 
                     #orf_list.append(match[a])
    return orf_list


############################################################################
#
# Function:     orfs_headers()
#
# Programmer:   Stephen Panossian
#
# Class:        BIFS 617 9040 Advanced Bioinformatics Summer 2020
#
# Assignment:   Group 3 Project
#
# Purpose:      Append frame number, genomic position, and ORF length values
#               to all header records from input sequence in FASTA format.
#
# Date:         July 28, 2020 / SP
#
# Version:      0.0
#               0.1 / 08-05-20 / SP modifed input variables and used regex
#                                   to find header
#
# Input:        headers  = header records from input FASTA file,
#                          list of strings
#               orf_list = list of 6 open reading frames
#
# Output:       headers  = header records appended with frame, position and
#                          length values, list of strings
#
# Notes:        FRAME_NO = global integer, frame number 1 through 6
#               POSITION = global integer, Forward 1, 2, 3/Reverse -1,-2,-3 
#               LEN_ORF  = global integer, number of bases in ORF
#
# References:
#
############################################################################

def orfs_headers(fasta1,orf_list):

    headers   = ["", "", "", "", "", ""] 
    pos_array = [1, 2, 3, -1, -2, -3]

    for line in fasta1:
        if (re.search(r"^\w>", line)):
            # found header record; append frame number, position, and length
            for i in range (1, 7):
                indx = i - 1
    
                FRAME_NO = i
    
                POSITION = pos_array[indx]

                LEN_ORF  = len(orf_list[indx])
    
                headers[indx] = line + " | FRAME = ", FRAME_NO + \
                        " POS = ", POSITION, + \
                        " LEN = ", LEN_ORF
            # end for
        # end if
    # end for

    return(headers)

############################################################################
#
# Function:     separate_orfs_into_codons()
#
# Programmer:   Ben Sparklin
#
# Class:        BIFS 617 9040 Advanced Bioinformatics Summer 2020
#
# Assignment:   Group Project
#
# Purpose:      Separate ORFs in orf_list into codons. Print ORFs in codon format with
#               no more than 15 codons per line.
#               
#               
#
# Date:         July 22, 2020 / BS
#               
#
# Version:      1.0
#               1.1/ 8-6-20 / BS Added a codon counter to ensure only 15 codons are printed per line.
#                   Each time the codon counter reaches 15, it resets to 0.
#               
#               
#
# Input:        orf_list
#
# Output:       codons
#
# References:   www.stackoverflow.com
#               www.python.org
#
############################################################################


def separate_orfs_into_codons():
    
    orf_list = determine_reading_frames_and_find_ORFs()
    
    for orf in orf_list: #for each ORF in the list of ORFs
        codons = "" #initilize codons variable
        codon_counter = 0 #initialize codon counter
        
        for i in range(0, len(orf), 3): #for each set of 3 nucleotides (1 codon) in the ORF
            codons += orf[i:i+3] + " " #add each codon to string with space between each codon
            codon_counter += 1 #each time a codon is added to string, add 1 to codon counter
            if codon_counter >= 15:
               codons += "\n" #add a newline to string each time 15 codons are recognized
               codon_counter = 0 #reset codon counter to 0
               
        print("ORF:\n", codons) #print the codons for each ORF with no more than 15 codons per line

################################################################################
#
# Program:          ORFS.py
#
# Programmers:      Mills, Mary (MM)
#                   Panossian, Stephen (SP)
#                   Sparklin, Ben (BS)
#                   Unnikumaran, Yekaterina (YU)
#
# Class:            BIFS 617 9040 Advanced Bioinformatics Summer 2020
#
# Assignment:       Group 3 Project
#
# Purpose:          To find all the open reading frames (ORFs) in an input
#                   DNA sequence in FASTA format.
#
# Date Submitted:   Aug. nn, 2020
#
# Date Due:         Aug. 11, 2020
#
# Version:          0.0 / SP (main)
#                   0.1 / SP - added prompt for minimum ORF length
#
# Notes:            None.
#
# References:       Jones, M. (2013) Chapter 8. Dictionaries.
#                   Python for Biologists.
#                   http://pythonforbiologists.com
#
#                   National Center for Biotechnology Information.
#                   https://www.ncbi.nlm.nih.gov
#
################################################################################

# Module(s) referenced
import re

def main(): 
    global headers
    global minimum_orf_length
    global fasta1
    global orf_list
    global accnum
    global seq
    
    headers = ["", "", "", "", "", ""] 

    # Introduce ORF program to user.
    print("\nWelcome to Group 3's Open Reading Frames (ORFS) Program!\n")
    print("ORFS will process five input DNA sequences and output all 6 " + \
          "open reading frames:" + \
          "\n\tFirst  nucleotide in codon, forward direction" + \
          "\n\tSecond nucleotide in codon, forward direction" + \
          "\n\tThird  nucleotide in codon, forward direction" + \
          "\n\tFirst  nucleotide in codon, reverse direction" + \
          "\n\tSecond nucleotide in codon, reverse direction" + \
          "\n\tThird  nucleotide in codon, reverse direction")

    
    #Collect input from user and parse for sequences and headers
    filesearchandread()
    
    # Read input file name
    read_input_file_name()
 
    # Read sequences from input file
    # (ensure that they are >= 50 characters)
    read_input_ORF_length()
    
    # Prompt user for minimum ORF length (>=50 bp)
    #while (minimum_orf_length < 50):
    #    str_min_ORF_len    = input("Enter your preferred ORF length in bp (minimum = 50)\n>>> ")
    #    minimum_orf_length = int(str_min_ORF_len)
    # end while

    # call computational function 1:
    #accnum, seq = filesearchandread()
    determine_reading_frames_and_find_ORFs()

    # Append headers list with ORF details.
    headers = orfs_headers(fasta1,orf_list)

    # Output ORFs to console
    separate_orfs_into_codons()

    print("\nGroup 3's Open Reading Frames (ORFS) Program has ended.\n")

main()
