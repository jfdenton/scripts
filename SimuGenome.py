#! /usr/bin/python
#This program will first analyze a genome, extracting nucleotide frequencies.  It will then create a new genome based upon the size of the input genome and the extracted frequencies.

#See readme for further notes - jfdenton

import cgi, random, sys, os, re

usage = "SimuGenome.py <input_file> <output_file>"

if len(sys.argv) != 3:  # require proper number of command line parameters
    print usage  #inform user of usage
else:
    input = sys.argv[1]
    output = sys.argv[2]

def genomeAnalyzer(input):
    bases = 0  #base counter
    genome = {'A':0,'G':0,'C':0,'T':0, 'N':0}
    A_, G_, C_, T_ = 0.0, 0.0, 0.0, 0.0  #init frequencies
    try:
        f = open(input)
        for line in f:
            for x in line:
                if x == 'A' or x == 'a':
                    genome['A'] += 1
                    bases += 1
                elif x == 'G' or x == 'g':
                    genome['G'] += 1
                    bases += 1
                elif x == 'C' or x == 'c':
                    genome['C'] += 1
                    bases += 1
                elif x == 'T' or x == 't':
                    genome['T'] += 1
                    bases += 1
                else:
                    genome['N'] += 0  # be wary of non-base characters and gaps
        print "Found %d total bases.." % bases 
        A_ = float(genome['A'])/float(bases)
        G_ = float(genome['G'])/float(bases)
        C_ = float(genome['C'])/float(bases)
        T_ = float(genome['T'])/float(bases)
        print "Nucleotide Composition:\nA: %.2f G: %.2f C: %.2f T:%.2f" % (A_, G_, C_, T_)
        Freq = A_ + G_ + C_ + T_  #sanity check, do our frequencies equal 1?
        if Freq != 1.0: 
            print "Frequencies do not add up, check your input file!"
        return A_, G_, C_, bases
        f.close()
    except:
        print "Failed to open input file!"  # check for typos in input file name

def genomeCreator(As, Gs, Cs, bases):
    o = open(output, 'a')  
    As = (As * 10000)    #consider 4 decimal places in random # comparison for accuracy
    Gs = (Gs * 10000) + As
    Cs = (Cs * 10000) + Gs
    for x in range(1,bases):   # Divide frequency into proportional spaces
        randy = random.randint(0,10000)  #Allow randy to place bases by proportion
        if randy <= As:
            o.write('A')
        elif randy <= Gs:
            o.write('G')
        elif randy <= Cs:
            o.write('C')
        else:
            o.write('T')
    print "Creating simulated genome in:", output
    o.close()

def diCounter(file, bases):
    gigantor = ''
    try:
        f = open(file)
    except:
        print "diCounter failed to open input files!"
        return
    eins = ['A', 'G', 'C', 'T']   # all nucleotides
    zwei = ['A', 'G', 'C', 'T']
    drei = {}   # a dictionary to house our frequencies
    for x in eins:   # loop through eins and zwei to create drei
        for y in zwei:
            drei[x+y] = 0.0  # frequency dictionary with all combinations
    for line in f:
        gigantor += line.strip()
    for x, y in drei.iteritems():
        regex = re.compile(str(x))
        matches = re.findall(regex, gigantor)
        drei[x] = "%.4f" % (float(len(matches))/float(bases))
    freq, count = '', 0
    for x,y in drei.iteritems():
        count += 1
        record = " %s: %.4f " % (x, float(y))
        freq += record
        if count % 2:
            freq += '\n'
    print freq

#----main program call-----#
print "Analyzing input genome..."
a, g, c, bases = genomeAnalyzer(input)
print "Determining di-nucleotide frequencies of real genome:"
diCounter(input, bases)
genomeCreator(a,g,c,bases)
print "Analyzing simulated genome.."
a,g,c, bases = genomeAnalyzer(output)
print "Determining di-nucleotide frequencies of simulated genome:"
diCounter(output, bases)
