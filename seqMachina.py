#! /usr/bin/python
import cgi, random, sys, os
from optparse import OptionParser

print "\n-=Welcome to Jimmy's Sequencing Machine Simulator=-"

usage = "Usage: python seqMachicha.py [options]\nUse -h or --help for further option information."

#---Main Variables & Containers---#

readLength = 100
reads = []
coords = {}
XCOV = 1.0

#---------------------------------#

def GenomeGrab():
    # This function extracts the sequence from the input file.
    gigantor = ''
    f = open(input_file)
    line = f.next()
    for line in f:
        gigantor += line.strip()
    m = len(gigantor)
    # returns length of input genome
    return m, gigantor

def Coverage( readNum ):                  # continuously calculate coverage until desired X
    cov = float(( ( float(readLength) * float(readNum)) / float(genomeLength) ))
    return cov

def Sequencer(m, g):              # accepts size of genome, and genome sequence string
    coverage, readNum, zero_bases = 0.0, 0, 0
    circular_ends = g[:201] + g[len(g)-200:-1]  # make string for reads spanning the ends
    while coverage <= XCOV:       # keep splicing reads until at chosen X coverage
        pos1 = random.randint(0,m)  # randomly choose splice site
        if pos1 < 100: read = circular_ends[pos1:pos1+100]  # if close to ends...
        elif pos1 > m-100: read = circular_ends[pos1:pos1-100] # use our spliced string
        else:                     # otherwise in middle of genome
            direction = random.randint(0,1)  # randomly choose direction (strand)
            if direction == 0:   # simulate minus strand
                pos2 = pos1-100
                read = g[pos1:pos2] 
                read = read[::-1]  # reverse the read
            else:
                pos2 = pos1+100
                read = g[pos1:pos2]
        if len(read) == 100:
            for x in range(pos1,pos2):    # record every position covered
                if x not in coords:
                    coords[x] = 0         # use index to keep track of unique occurences
            reads.append(read)
            readNum += 1
            coverage = Coverage(readNum)
    uncovered =  (float(genomeLength) - float(len(coords)) )/ float(genomeLength) 
    print 'Total Positions: %d \nProportion Not Covered: %.4f' % (len(coords),uncovered )

def Output():
    o = open(output_file, 'a')            # simple output function
    i = 1
    for x in reads:
        record = '>read%s\n%s\n' % (i,x)  # fasta-like format
        i += 1
        o.write(record)
    print 'Output created in:', output_file

#---------Main---------#

    #----Options----#
parser = OptionParser(usage=usage)
parser.add_option("-i", "--INPUT", dest="input_file", help="The path to the input sequence file.")
parser.add_option("-o", "--OUTPUT", dest="output_file", help="The name of the desired output file for the sequence alignment.")
parser.add_option("-c", "--COVERAGE", dest="coverage", help="The average coverage of the genome you'd like to achive.")

(options, args) = parser.parse_args()

if (len(sys.argv) < 4):
    print(usage)
    sys.exit()

if (not os.path.exists(options.input_file)):
    print("Could not open input file!" % ( options.input_file))
    sys.exit()

input_file = options.input_file
output_file = options.output_file
if (options.coverage):
    XCOV = float(options.coverage)

    #----Program Calls----#

genomeLength, genome = GenomeGrab()
#grab genome sequence and length

Sequencer(genomeLength, genome)
#simulated sequencing machine will splice reads out and compute uncovered frequency

Output()

    #----------------------#



