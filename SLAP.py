#! /usr/bin/python
import cgi, random, sys, os
from optparse import OptionParser

print "\nWelcome to Jimmy's Smith-Waterman Local Alignment Program in Python [SLAP]"

usage = "Usage: python SLAP.py [options]\nUse -h or --help for further option information."

#---Main Variables & Containers---#

seqs, matrix, seqn = ['',''], [], 0   # useful global variables
M,N = 0,0                             # lengths of our two sequences
gap_open, gap_extend = -11, -1        # penalties
dynosaur, delete, insert = [], [], [] # dynamic programming arrays
tracetrail = []                       # tracetrail - path saving array

# Create protein index to relate to blossom matrix
PI = {'A':0, 'R':1, 'N':2, 'D':3, 'C':4, 'Q':5, 'E':6, 'G':7, 'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16, 'W':17, 'Y':18, 'V':19, 'B':20, 'J':21, 'Z':22, 'X':23, '*':24}


#---------------------------------#

def SeqGrab():
    # This function extracts the sequences from the input file.
    f = open(input_file)
    gigantor = ''
    for line in f:
        if line.startswith('>'):
            gigantor += '-'
        else:
            gigantor += line.strip()
    gigantor = gigantor.split('-')
    seqs[0], seqs[1] = gigantor[1], gigantor[2]   
    m,n = len(seqs[seqn]), len(seqs[seqn+1])
    # returns length of sequences to global variables M,N
    return m, n

def MatrixGen():
    # This function generates protein scoring matrix from BLOSUM62 file.
    try:
        blosum = open('BLOSUM62')
        for line in blosum:   #fill up the matrix with appropriate values
            row = []
            try:
                test = int(line[3])  #test to look for integers
                line = line[1:].strip()
                splitter = line.split(' ')
                for x in splitter:
                    if x != '':      #watch for double spaces
                        row.append(x)
                matrix.append(row)
            except:        
                continue  # essentially skip the first line     
    except:
        print "Could not open Blossum Matrix File ('BLOSUM62')!"
        sys.exit()

def ArrayBuilder(array):
    #simple function to generate an array of arrays and initliaize to 0
    for x in range(M+1):
        y = [0 for i in range(N+1)]
        array.append(y)
        
def DynamicBuilder():
    # Main algorithm scores alignment, making use of 3 arrays for comparison, and 1 array for path saving
    m,n = M+1,N+1
    Max, iMax, jMax = 0,0,0   # Initialize max values to 0
    for i in range(1,m):      # Now use dynamic programming to recursively fill the arrays
        for j in range(1,n):
            delete[i][j] = max(dynosaur[i-1][j]+GapPen(seqn+1), delete[i-1][j]+gap_extend)
            insert[i][j] = max(dynosaur[i][j-1]+GapPen(seqn+1), insert[i][j-1]+gap_extend)
            blosum = matrix[ PI[seqs[seqn][i-1]] ][ PI[seqs[seqn+1][j-1]] ]  # get match score from BLOSUM62
            match, gap1, gap2= ( dynosaur[i-1][j-1] + int(blosum) ), delete[i][j], insert[i][j]
            dynosaur[i][j] = max(0, gap1, gap2, match)             # Choose max, magic 0 prevents negatives
            if dynosaur[i][j] == seqn: tracetrail[i][j] = seqn       # also save markers for traceback
            if dynosaur[i][j] == gap1: tracetrail[i][j] = seqn+1     # deletion
            if dynosaur[i][j] == gap2: tracetrail[i][j] = seqn+2     # insertion
            if dynosaur[i][j] == match: tracetrail[i][j] = seqn+3    # match / mismatch
            if dynosaur[i][j] >= Max: iMax, jMax, Max = i, j, dynosaur[i][j]
    return int(iMax), int(jMax)                              # return max of i,j so we can trace our path

def Traceback(i,j):                                          # as we scored, we saved our path to tracetrail
    seqAlign1, seqAlign2 = '',''
    while tracetrail[i][j] != 0:                             # while not at end of path, count back
        if tracetrail[i][j] == 1:
            seqAlign1 = seqAlign1 + seqs[seqn][i-1]          # saving alignment as we go
            seqAlign2 = seqAlign2 + '-'
            i -= 1
        elif tracetrail[i][j] == 2:                          # follow our markers backwards
            seqAlign1 = seqAlign1 + '-'
            seqAlign2 += seqs[seqn+1][j-1]
            j -= 1                                     
        elif tracetrail[i][j] == 3:
            seqAlign1 = seqAlign1 + seqs[seqn][i-1]
            seqAlign2 = seqAlign2 + seqs[seqn+1][j-1]
            i -= 1
            j -= 1
    return seqAlign1[::-1], seqAlign2[::-1]                  # flip aligns as we return

def GapPen(length):
    gap_penalty = (gap_open + length * gap_extend)
    return gap_penalty

def Scoring(x,y):                                             # accurate scoring, create hit_line for reporting
    hit_line, score, gapStart = '', 0, False                  
    for i in range(len(max(x,y))):                            # for the length of the alignment
        if x[i] == y[i]:
            hit_line += '|'                                   # | will connect matches
            score += int(matrix[ PI[x[i]] ][ PI[y[i]] ])      # get BLOSUM62 match score
            gapStart = False
        elif x[i] != y[i] and x[i] != '-' and y[i] != '-':
            hit_line += '.'                                   # use . for mismatch
            score += int(matrix[ PI[x[i]] ][ PI[y[i]] ])      # calculate mismatch BLOSUM62 score
            gapStart = False
        elif x[i] == '-' or y[i] == '-':
            hit_line += ' '
            if gapStart == False:                             # affine gap penalties, -11 open, -1 extend
                score += gap_open
                gapStart = True
            elif gapStart == True:
                score += gap_extend
    return hit_line, score

def Outputter(align1, align2, hit_line, score):
    o = open(output_file, 'a')                                # finally, report results to output file
    record = 'seq%d:\t%s\n\t%s\nseq%d:\t%s\n' % (seqn+1, align1,hit_line,seqn+2, align2)
    print 'Alignment Score:',score                            # report alignment score also
    print 'Output Created in:',output_file
    if file_format == 'FASTA':
        record = ">seq1\n" + align1 + '\n' + ">seq2\n" + align2 + '\n'
        o.write(record)
    else:
        o.write(record)
        

#---------Main---------#

    #----Options----#
parser = OptionParser(usage=usage)
parser.add_option("-i", "--INPUT", dest="input_file", help="The path to the input sequence file.")
parser.add_option("-f", "--FORMAT", dest="file_format", help="The desired output format of the alignment.[DEFAULT | FASTA]")
parser.add_option("-o", "--OUTPUT", dest="output_file", help="The name of the desired output file for the sequence alignment.")

(options, args) = parser.parse_args()

if (len(sys.argv) < 5):
    print(usage)
    sys.exit()

if (not os.path.exists(options.input_file)):
    print("Could not open input file!" % ( options.input_file))
    sys.exit()

input_file = options.input_file
output_file = options.output_file
file_format = options.file_format

    #----Program Calls----#
    
# extract sequences from input fasta, store, and return their length
M,N = SeqGrab()
# construct the matrix from BLOSUM62 input file
MatrixGen()

# Establish appropriately sized dynamic programming arrays. Initialize values to 0
ArrayBuilder(tracetrail)
ArrayBuilder(dynosaur)
ArrayBuilder(delete)
ArrayBuilder(insert)

# Fill in DP arrays and Traceback table, fill in recursively, return final index positions.
i,j = DynamicBuilder()

# Traceback through tracetrail array, using our 4 markers, building the sequence as we go, flip them and return
align1, align2 = Traceback(i,j)

# Toss alignments through more accurate scoring function, and create hit_line for default output format
hitLine, score = Scoring(align1,align2)

# Finally Output based on file_format option, print alignment score
Outputter(align1, align2, hitLine, score)

    #----------------------#



