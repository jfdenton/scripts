#! /usr/bin/python
#This robust program accepts as input a protein sequence in fasta format, and will compute the six possible protein configurations, three shifts for forward and three for reverse.  Finally it will report these proteins in a somewhat fasta like file.
# See readme for more notes ~ jfdenton
import cgi, random, sys, os

usage = "BioTranslator.py <input> <output>\n"

AAs, Base1, Base2, Base3 = '', '', '', ''
code, codonated, translated, gigantor = {}, {}, {}, []

def bioTranslate():
    entries, shift = 0,0
    try:
        f = open(input)
        o = open(output, 'a')
    except:
        print "Failed to open file!"  # typed file names correctly check
        return
    for line in f:
        if line[0] == '>':
            entries += 1      # keep track of number of entries
            name = line[:].strip()+"| protein list ID:"  # create unique fasta id
            line = f.next()
            line = line.replace('N', 'G')  # not handling well input with N nucleotides
            gigantor.append([name+' '+str(entries),line.strip()])  # add number for sorting later
        else:   # lets make sure we dont miss sequences that extend beyond 1 line, unlike sampleseq.txt
            if line[0] != '>' and line[0] != ' ' and line[0] != '':
                line = line.replace('N','T') # again, just in case input has Ns, convert
                gigantor[entries-1]+=line
    for x in gigantor:    # forward strand codon extraction
        codonated[x[0]] = [[],[],[],[],[],[]]  #establish dictionary entry, 2 strands * 3 shifts = 6 slots
        for y in range(3):
            codonated[x[0]][shift] = codonExtract(x[1], shift)  # extract codons for all 3 shifts
            shift += 1
        shift = 0
    shift = 0
    for x in gigantor:    # reverse strand codon extraction
        translated[x[0]] = [[],[],[],[],[],[]] # looking ahead, need space for the translated sequence
        for y in range(3):   # update diciontary entries 4-6
            codonated[x[0]][shift+3] = codonExtract(reverseComp(x[1]), shift)
            shift += 1     # note the call above to reverseComp within the call to codonExtract
        shift = 0
    shift = 0
    for x, y in codonated.iteritems():  #iterate through dictionary of entries, codons
        for items in y:
            sequence = ''
            for triplets in items:
                protein = code[triplets]  #toss codons through our protein code dictionary
                sequence += protein
            translated[x][shift]=sequence  #update translated dictionary for all 6 sequences
            shift += 1
        shift = 0
    final = translated.items()  # sorting dictionaries is never fun, convert to a list
    final.sort(key=lambda translated:translated[1], reverse=True)
    final.sort( key=lambda translated:translated[0] )
    for x,y in final:   #Now finally iterate through our translated proteins and report
        record = ">%s\n>Forward\n%s\n%s\n%s\n>Reverse\n%s\n%s\n%s\n" % (x,y[0],y[1],y[2],y[3],y[4],y[5])
        o.write(record)      #output file is somewhat fasta like with names, forward, and reverse listed
    print "Found %d entries, computed 6 protein sequences for each, in file: %s" % (entries, output)
    f.close()
    o.close()
            
           
def codonExtract(x, shift):   #iterates over sequence gathering codons
    codons = []                 #accepts shift as well as sequence
    start, end  = 0, 3
    bases = len(x)
    for y in range(0, bases/3):
        codon = x[start+shift:end+shift]
        start += 3
        end += 3
        if len(codon) == 3:
            codons.append(codon)
    start, end = 0, 3
    return codons              # return codons to be added to codonated dictionary

def reverseComp(seq):   # Returns reverse complement of sequence
   rev_seq = ''
   for nucleotide in seq:  # switch those...
       if nucleotide == 'T':
           rev_seq += 'A'
       if nucleotide == 'G':
           rev_seq += 'C'
       if nucleotide == 'A':
           rev_seq += 'T'
       if nucleotide == 'C':
           rev_seq += 'G'
   rev_seq = rev_seq[::-1]  # reverse it...
   return rev_seq

def codonBuilder():   # Attempt to parse codon.txt, if it is present
    try:
        c = open('codon.txt')
        for line in c:
            if 'AAs' in line:
                AAs = line.strip()
                split = AAs.split('= ')
                AAs = split[1]
            if 'Base1' in line:
                Base1 = line.strip()
                split = Base1.split('= ')
                Base1 = split[1]
            if 'Base2' in line:
                Base2 = line.strip()
                split = Base2.split('= ')
                Base2 = split[1]
            if 'Base3' in line:
                Base3 = line.strip()
                split = Base3.split('= ')
                Base3 = split[1]
    except:           # If its not present, manually create the lists
        print "Failed to open codon.txt, manually creating the lists!"
        AAs = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        Base1='TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
        Base2='TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
        Base3='TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
    for x in range(0,len(AAs)):   # create reference dictionary for codon -> protein conversion
        code[Base1[x]+Base2[x]+Base3[x]] = AAs[x]
        
    
#----Main Program Calls---#

if len(sys.argv) != 3:   # check for proper command line arguments
    print usage          # remind user of usage
else:
    input = str(sys.argv[1])
    output = str(sys.argv[2])                  
    codonBuilder()      # build codon table from parsing codon.txt, or from strings
    bioTranslate()      # main function call

