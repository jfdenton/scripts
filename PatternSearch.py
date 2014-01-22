#! /usr/bin/python
#This is a basic version of my PatternSeek (PatternSearch2.py). It only reports occurences, while the updated and perhaps overzealous version makes an effort to report the line they were found upon as well as the positions within the line.  Both scripts accept standard, python, or Prosite format for pattern input. 
#~ jfdenton

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
    patternMatch(regex)
    
#Main Function: For a given input file, this function will seek to match pattern input by the user.  It will report the occurences, and the locations of the occurences.  

def patternMatch(regex):
    gigantor, matches = '', []
    just_checking, blahbu = [], []
    f = open(input)
    for line in f:
        gigantor += line.strip()   #build a string from the lines, remove newlines
        match = re.findall(regex, line)
        if len(match) > 0:
            matches.append(match[0].strip())
            match = []
    just_checking = re.findall(regex, gigantor)  #make sure we didnt miss any between lines
    for x in just_checking: 
        if x not in matches:   #if we did, keep track of that count too
            blahbu.append(x)
    
    print "Jimmy's patternMatch found %d matches!" % (len(matches)+len(blahbu))
    f.close()

#------Main Program Calls-------#

if pattern[0] == '(' or pattern[-1] == ')':
    pattern = pattern[1:-1]
if 'x(' in pattern or '-' in pattern: #Prosite format Detected!
    patternToRE(pattern) #Call the converter
else:
    regex = re.compile(pattern)  #insure regex format
    patternMatch(regex)          #call main function
