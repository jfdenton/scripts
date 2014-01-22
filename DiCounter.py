#! /usr/bin/python
#This is a standalone version of diCounter function within
#SimuGenom.py - It accepts as input two files, the real and simulated raw genome sequence files. It will report frequency of all possible di-nucleotides.

import cgi, random, sys, os, re

usage = "DiCounter.py <input_file:real> <input_file:simulated>"

def diCounter(gigantor, bases):
    eins = ['A', 'G', 'C', 'T']   # all nucleotides
    zwei = ['A', 'G', 'C', 'T']
    drei = {}   # a dictionary to house our frequencies
    for x in eins:   # loop through eins and zwei to create drei
        for y in zwei:
            drei[x+y] = 0.0  # frequency dictionary with all combinations
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

#--- Main Program Calls --- #


if len(sys.argv) != 3:  # require proper number of command line parameters
    print usage  #inform user of usage
else:
    input = sys.argv[1]
    output = sys.argv[2]
    try:
        file = open(input)
    except:
        print "Couldn't open real genome file!"
    bigantor = ''
    for line in file:
        bigantor += line.strip()
    diCounter(bigantor, len(bigantor))
    print "-----------------------------------"
    bigantor = ''
    try:
        file = open(output)
    except:
        print "Couldn't open simulated genome file!"
    for line in file:
        bigantor += line.strip()
    diCounter(bigantor, len(bigantor))
