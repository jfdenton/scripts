#! /usr/bin/python

import cgi, random, sys, os, re

usage = "PatternSearch.py <pattern_to_match> <input_file>"

if len(sys.argv) != 3:  # require command line parameters
    print usage  # inform user of proper command line usage
else:
    pattern = sys.argv[1]  # assign inputs to pattern and I/O
    input = sys.argv[2]

#I've included this function to allow as pattern input ProSite format

def patternToRE(pattern):
    pattern = pattern.replace(' ','')
    pattern = pattern.replace('-','')
    pattern = pattern.replace('x','.')
    pattern = pattern.replace('{','[^')
    pattern = pattern.replace('}',']')
    pattern = pattern.replace('(','{')
    pattern = pattern.replace(')','}')
    pattern = pattern.replace('>','$')
    regex = re.compile(pattern)
    patternSeek(regex)
    
#Main Function: For a given input file, this function will seek to match pattern input by the user.  It will report the occurences, and the locations of the occurences.  

def patternSeek(regex):
    maybe, matched, lineCount, matches, gigantor = False, False, 0, [], ''
    record = True
    f = open(input)
    for line in f:
        gigantor+=line.strip()
        lineCount += 1   #how I report which line the match was on
        #if lineCount == number: print line  #if you wish to test the coordinates
        match = re.findall(regex, line)  #change the number var to line you wish to see
        if len(match) > 0:   #if we find a match, start iterating through the line
            pos = 0
            coords = [0]
            for x in line:  #to extract correct coordinates, 2 state variables are used
                pos += 1
                if x == match[0][0] and maybe == False:
                    coords[0] = pos
                    maybe = True
                if maybe == True:
                    if x == match[0][1]:
                        matched = True
                    else:
                        maybe = False
                if matched == True:
                    if x == match[0][-1]:
                        record = True
            if record == True:  #if 3 proper positions are found, record the match coords
                position = (coords[0], coords[0]+7)
                matches.append((match,position,lineCount))
                coords, match, position = [], [], ()  # reset stats
                maybe, matched = False, False
    matchez = re.findall(regex, gigantor)

    #Report complete pattern match findings
    split = input.split('.')
    prefix = split[0]
    output = prefix+"_vs_"+pattern+".txt"  #create output file for coordinates
    o = open(output, 'a')
    for x in matches:
        record = "%s on line %d at position (%d,%d).\n" % (x[0][0],x[2],x[1][0],x[1][1])
        o.write(record)
    print "\nJimmy's patternSeek found a total of %d matches!" % (len(matchez))
    print "Results were placed in %s in the current directory." % (output)
    f.close()
    o.close

#------Main Program Calls-------#


if pattern[0] == '(' and pattern[-1] == ')':  #check for parantheses in pattern input
    pattern = pattern[1:-1]   #bash will likely scream if you try them anyway

if 'x(' in pattern or '-' in pattern: #Prosite format Detected!
    patternToRE(pattern) #Call the converter
else:
    regex = re.compile(pattern)  #insure regex format
    patternSeek(regex)          #call main function


