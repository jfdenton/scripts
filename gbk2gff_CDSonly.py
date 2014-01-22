#! /usr/bin/python

import sys, os, re
from optparse import OptionParser

usage = "Usage: python gbk_parser.py [options]\nUse -h or --help for further option information."

#EXAMPLE1 USAGE: python gbk_parser.py -i ame_ref_Amel_4.5_chrLG1.gbk
#EXAMPLE2 USAGE: python gbk_parsery.py -i ame_ref_Amel_4.5_chrLG1.gbk -o Amel_4.5_chrLG1 -c 1 -v 1

states = [ "LOCUS" , "DEFINITION", "ACCESSION" , "VERSION" , "DBLINK" , "KEYWORDS", "SOURCE", "ORGANISM", "REFERENCE", "COMMENT" , "FEATURES" , "ORIGIN" ]

def convertGBKtoGFF3(file_name):
    # This function will read over the input genbank file and create separate gff3 files for each locus
    dbx, multi, multi_seq, feat_start = False, False, False, False
    featNum, state, sub_state = 0, "LOCUS", "source_info"
    loci, features, ids, converted = {}, {}, {}, []
    gff_header, comment, note = "##gff-version 3\n", ";comment1=", ""
    f = open(file_name)
    for line in f:
        splitor = line.split()
        if line.startswith("//"):    #end of entry, write out the locus features and re-init
            multi, dbx = False, False
            features[0] = gff_header
            loci[locus] = features
            o = open(locus + '.gff3', 'w')  # user could change this output prefix (GFFOUT)
            o.write(gff_header)
            for x in range(1,len(loci[locus])):
                if x in loci[locus]:
                    Ar = loci[locus][x]
                    record = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Ar[0], Ar[1], Ar[2], Ar[3], Ar[4], Ar[5], Ar[6], Ar[7], Ar[8])
                    o.write(record)
            o.close()
            converted.append(locus + '.gff3')
            features = {}
            gff_header, comment, note = "##gff-version 3\n", ";comment1", ""
            gff_fields = ["locus", "Genbank", "feature", "start", "stop", ".", "strand", "frame", "attr"]
            feat_start = False
            featNum = 1

        #state update
        for each in states:
            if line.startswith(each) or splitor[0] == each:
                state = each
        if state == "FEATURES": feat_start = True
        if state == "ORIGIN": feat_start = False
        
        #state instructions
        if state == "LOCUS":
            splitee = re.split(r'\s{2,}',line)
            locus = splitee[1]
            splitee = line.split()
            date = splitee[7]
        elif state == "DEFINITION":      
            note += ' ' + line[12:].strip()    
        elif state == "COMMENT":
            comment += ' ' + line[12:].strip()
        elif state == "FEATURES" and feat_start:
            # update sub_state
            if feat_start and not line.startswith("FEATURES"):
                # booleans for grabbing multi-line features
                if not re.match(r'^(\s){21}/db_x', line): dbx = False
                if multi_seq == True and re.match(r'^(\s){21}/', line): multi_seq = False
                if multi == True and re.match(r'^(\s){21}/', line): multi = False

            # construct gff header and region entry
                if re.match(r'^(\s){5}source', line):
                    featNum += 1
                    sub_state = "source"
                    splitee = line[21:].strip().split('..')
                    features[featNum] = [locus, "Genbank", "region", "start", "stop", ".", "+", ".", "attr"]
                    features[featNum][3] = splitee[0]
                    features[featNum][4] = splitee[1]
                    gff_header += "# sequence-region " + locus + ' ' + splitee[0] + ' ' + splitee[1] + '\n'
                    features[featNum][8] = "ID=" + locus + ';Name=' + locus
            # gene state init
                elif re.match(r'^(\s){5}gene', line):
                    featNum += 1
                    sub_state = "gene"
                    splitee = line[21:].strip().split('..')
                    if splitee[0].startswith('complement'):
                        splitoo = splitee[0].split('(')
                        splitee[0] = splitoo[1]
                        splitee[1] = splitee[1][:-1]
                        strand = "-"
                    else: strand = "+"
                    features[featNum] = [locus, "Genbank", "gene", "start", "stop", ".", strand, ".", "attr"]
                    features[featNum][3] = splitee[0]
                    features[featNum][4] = splitee[1].strip()
                    features[featNum][8] = "ID="
            # mRNA state init
                elif re.match(r'^(\s){5}mRNA', line):
                    sub_state = "mRNA"
                    if line[21:].startswith('complement(join('):
                        strand = '-'
                        coords = line[21:].strip().split('join(')[1]
                        coords = coords.split(')')[0]
                        multi_seq = True
                    elif line[21:].startswith('join('):
                        strand = '+'
                        splitor = line[21:].strip().split('join(')[1]
                        coords = splitor.split(')')[0]
                        multi_seq = True
                    elif line[21:].startswith('complement('):
                        strand = '-'
                        coords = line[21:].strip().split('complement(')[1]
                        coords = coords.split(')')[0]
                    else:
                        strand = "+"
                        coords = line[21:].strip()
                    featNum += 1
                    attr = ""
                    features[featNum] = [locus, "Genbank", "mRNA", "start", "stop", ".", strand, ".", attr]
            # CDS state init
                elif re.match(r'^(\s){5}CDS', line):
                    sub_state = "CDS"
                    if line[21:].startswith('complement(join('):
                        strand = '-'
                        coords = line[21:].strip().split('join(')[1]
                        coords = coords.split(')')[0]
                        multi_seq = True
                    elif line[21:].startswith('join('):
                        strand = '+'
                        splitor = line[21:].strip().split('join(')[1]
                        coords = splitor.split(')')[0]
                        multi_seq = True
                    elif line[21:].startswith('complement('):
                        strand = '-'
                        coords = line[21:].strip().split('complement(')[1]
                        coords = coords.split(')')[0]
                    else:
                        strand = "+"
                        coords = line[21:].strip()
            # Gap state init
                elif re.match(r'^(\s){5}gap', line):
                    sub_state = "gap"
                    featNum += 1
                    splitee = line[21:].strip().split('..')
                    if splitee[0].startswith('complement'):
                        splitoo = splitee[0].split('(')
                        splitee[0] = splitoo[1]
                        splitee[1] = splitee[1][:-1]
                        strand = "-"
                    else: strand = "+"
                    features[featNum] = [locus, "Genbank", "gap", "start", "stop", ".", strand, ".", "attr"]
                    features[featNum][3] = splitee[0]
                    features[featNum][4] = splitee[1].strip()
                    features[featNum][8] = "ID=Genbank:gap:" + locus + ":" + splitee[0] + ":" + splitee[1] + ';'
                    features[featNum][8] += "Name=Genbank:gap:" + locus + ":" + splitee[0] + ":" + splitee[1] + ';'
            # miscRNA state init
                elif re.match(r'^(\s){5}misc_RNA', line):
                    sub_state = "misc_RNA"
            # tRNA - not currently implemented, later will add user option to include nc features
                elif re.match(r'^(\s){5}tRNA', line):
                    sub_state = "tRNA"
            # ncRNA - not currently implemented, later will add user option to include nc features
                elif re.match(r'^(\s){5}ncRNA', line):
                    sub_state = "ncRNA"
            # for debugging sub_state feature headers we did not consider
                else:
                    if re.match(r'^(\s){5}(\w+)(\s+)', line):
                        sub_state = line[:21].strip()
                        # We are only concerned with a subset of features, but should eventually handle any and all features, but allow user to choose (nc/c, etc)
                        # full_feats_list = ["source", "gene", "mRNA", "CDS", "gap", "misc_RNA", "ORIGIN", "tRNA", "Taxon", "protein_id", "GI", "translation", "complement", "promoter", "TATA_signal", "5'UTR", "3'UTR", "protein_bind", "exon", "repeat_region", "misc_feature", "rep_origin", "RBS",  "intron", "polyA_signal", "V_region", "sig_peptide", "J_segment", "C_region", "mat_peptide", "-35_signal" ]
                         
            # region sub_state instructions
            if sub_state == "source" and multi == False:
                if re.match(r'^\s{21}\/organism', line):
                    splitee = line.split('"')
                    org = splitee[1]
                    #splitoo = re.split('\W+', splitee[1])
                    gff_header += '# organism ' + org + '\n'
                    gff_header += '# date ' + date + '\n'
                    gff_header += '# Note' + note + '\n'
                    features[featNum][8] += ";organism=" + org
                elif re.match(r'^\s{21}\/db_xref=', line):
                    splitee = line.split('/db_xref="')
                    splitoo = splitee[1].split('"')
                    features[featNum][8] +=  ';Dbxref=' +splitoo[0] + ';Note=' + note
                    features[featNum][8] += comment + ';date=' + date
                elif re.match(r'^\s{21}\/note=', line):
                    multi = True
                elif re.match(r'^\s{21}\/mol_type', line):
                    splitee = line.split('"')
                    features[featNum][8] += ";mol_type=" + splitee[1]
                elif re.match(r'^\s{21}\/strain', line):
                    splitee = line.split('"')
                    features[featNum][8] += ";strain=" + splitee[1]

            # GENE sub_state instructions                   
            elif sub_state == "gene" and multi == False:
                if re.match(r'^\s{21}\/gene', line):
                    splitee = line.split('"')
                    ids[splitee[1]] = {'r':0, 't':0, 'p':0}
                    #print splitee[1]
                    features[featNum][8] += splitee[1] + ';Name=' + splitee[1] + ';gene=' + splitee[1]
                elif re.match(r'^\s{21}\/note', line):
                    multi = True
                    splitee = line[21:].strip().split('"')
                    features[featNum][8] += ';Note=' + splitee[1] 
                elif re.match(r'^\s{21}\/db_xref', line) and dbx == False:
                    features[featNum][8] += ';Dbxref=' + line[21:].split('"')[1]
                    dbx = True
                elif re.match(r'^\s{21}\/db_xref', line) and dbx == True:
                    features[featNum][8] += ',' + line[21:].split('"')[1]
            elif sub_state == "gene" and multi == True:
                if not re.match(r'^\s{21}\/', line): features[featNum][8] += line[21:].strip()
                else: multi = False

            # GAP sub_state        
            elif sub_state == "gap":
                if re.match(r'^\s{21}\/estimated_length', line):
                    splitee = line.strip().split('/')
                    features[featNum][8] += splitee[1]

            # mRNA sub_state
                    
            elif sub_state == 'mRNA' and multi_seq == False:
                if re.match(r'^\s{21}\/gene', line):
                    splitee = line.split('"')
                    gene = splitee[1]
                    ids[gene]['t']+=1  #gene id update
                    features[featNum][8] = 'ID=' + gene + '.t'+ str(ids[gene]['t']) + ';Name='+gene+'.t' + str(ids[gene]['t'])
                    while ')' in coords: coords = coords.split(')')[0]
                    splitor = coords.split(',')
                    maxy, miny = 0, 1000000000
                    mFeat = featNum
                    attr = "Parent=" + gene + '.t'+ str(ids[gene]['t']) + ';gene='+gene #parent id creation
                    for x in splitor:  # make exon entries
                        featNum += 1
                        #Catches some weird carrots in some gbk files inserted with start, stop coords, a regex would be more concise here (\d+)(?:...)(\d+)
                        if '<' in x: x = x.split('<')[1]
                        elif '>' in x: x = x.split('>')[0] + x.split('>')[1]
                        if '<' in x: x = x.split('<')[1]
                        elif '>' in x: x = x.split('>')[0] + x.split('>')[1]
                        start, stop = x.split('..')[0], x.split('..')[1]
                        frame = (int(start) - int(stop)) % 3
                        miny = min(min(int(start), int(stop)), miny)
                        maxy = max(max(int(start), int(stop)), maxy)
                        features[featNum] = [locus, "Genbank", "exon", start, stop, ".", strand, frame, attr]
                    features[mFeat][3], features[mFeat][4] = miny, maxy
                elif re.match(r'^\s{21}\/note', line):
                    multi = True
                    splitee = line[21:].strip().split('"')
                    features[featNum][8] += ';Note=' + splitee[1] 
                elif re.match(r'^\s{21}\/db_xref', line) and dbx == False:
                    features[featNum][8] += ';Dbxref=' + line[21:].split('"')[1]
                    dbx = True
                elif re.match(r'^\s{21}\/db_xref', line) and dbx == True:
                    features[featNum][8] += ',' + line[21:].split('"')[1]
            elif sub_state == "mRNA" and multi == True and multi_seq == False:
                if not re.match(r'^\s{21}\/', line): features[featNum][8] += line[21:].strip()
                else: multi = False

            elif sub_state == "mRNA" and multi_seq == True and multi == False and '(' not in line:
                coords += line[21:].strip()

            # CDS sub_state instructions
            elif sub_state == 'CDS' and multi_seq == False:
                
            elif sub_state == 'CDS' and multi_seq == False and multi == True:
            
            elif sub_state == 'CDS' and multi_seq == True and multi == False and '(' not in line::

            # miscRNA sub_state instructions         

    return converted
    
                    
def validateGFF3(x):
    never = 1 + 1
    #validation subroutines here

                
#----Options----#

parser = OptionParser(usage=usage)
convert = 0
validate = 0
parser.add_option("-i", "--input", dest="input_file", help="The path to the input genbank or gff3 file.")
parser.add_option("-o", "--output", dest="gff_out", help="The path to the output gff3 file prefix.")
parser.add_option("-c", "--convert", dest="convert", help="Convert genbank to gff3? -c 1 [default -c 0]")
parser.add_option("-v", "--depth", dest="validate", help="Run validation function? -v 1 [default -v 0]")
(options, args) = parser.parse_args()

if (len(sys.argv) < 2):
    print(usage)
    sys.exit()

if (not os.path.exists(options.input_file)):
    print("Could not open input file!" % ( options.input_file))
    sys.exit()

input_file = options.input_file

#--main function calls--#

#check option parameters before calling conversion script, or check if input was gbk
if convert or input_file.split('.')[-1] == 'gbk':
    outFiles = convertGBKtoGFF3(input_file)
    print "Converting GBK %s generated %d files:\n" % (input_file, len(outFiles))
    print outFiles

if validate:
    validateGFF3(x)


#This can be modified to test values held in memory from the genbank file, for debugging
query = """
queries = raw_input("Would you like to query this gbk file?")
if queries == "y" or queries == "yes" or queries == "Y" or queries == "Yes":
    print loci.keys()
    while query != "bye":
        query = raw_input("Enter the loci to query or type 'bye' to exit: ")
        if query in loci.keys():
            query = raw_input("Enter the feature to query or type 'bye' to exit: ")
        if query.upper() in states:
            print states[str(query.upper())].keys()
            source_pick = raw_input("Choose a source field :")
            if str(source_pick) in states[query.upper()]:
                print states[query.upper()][source_pick]
            else: print "I don't seem to contain that field.\n"

print "Thank you, come again."
  
"""

