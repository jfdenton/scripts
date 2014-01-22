#! /usr/bin/python

import cgi, random, sys, os, math
import numpy as np

print "\n-=Welcome to Jimmy's Coverage Calculator v2.0=-"
usage = "\nUsage: coverage <input_genes> <input_pileup> <output_genes>\n"

if len(sys.argv) < 3:
    print usage
    exit()
else:
    gene_in = str(sys.argv[1])
    pileup_file = str(sys.argv[2])
    gene_out = str(sys.argv[3])
    cov_out = str(sys.argv[3])+'cov.out'
    covBins = str(sys.argv[3])+'.bins'
    binSize = 1

gene_dict = {}
cov_dict = {}
exons = {}

def dict_gen():
    for x in range(19000):
        scaff = 'scaffold:'+str(x)
        gene_dict[scaff] = []
    f = open(gene_in)
    for line in f:
        splitter = line.strip().split(' ')
        geneNum = splitter[0]
        scaffSplit = splitter[1].split(':')
        #scaff = scaffSplit[1]
        scaff = str(splitter[1])
        start = int(splitter[2])
        end = int(splitter[3])
        #if scaff not in gene_dict:
        #    gene_dict[scaff] = []
        gene_dict[scaff].append( [start,end,geneNum] )

def coverage():
    f = open(pileup_file)
    o = open(gene_out, 'w')
    c = open(cov_out, 'w')
    start = False
    for line in f:
        splitter = line.split('\t')
        chr = str(splitter[0])
        if start == False:
            start = True
            record = "%s\n" % chr
            o.write(record)
            #print 'Chromosome %s has %d exons to check' % (chr, ) )
        splitter = line.split('\t')
        chr = str(splitter[0])
        pos = int(splitter[1])
        coverage = int(splitter[3])
        for exon in gene_dict[chr]:
            if pos >= exon[0]  and pos <= exon[1] or pos >= exon[1] or pos <= exon[0]:
                record = "%d %d\n" % (exon[0], exon[1])
                if exon in gene_dict[chr]: 
                    gene_dict[chr].remove(exon)
                    o.write(record)
                gene_name = chr +'_'+ exon[2]
                cov_dict[gene_name] = cov_dict.get(gene_name,0.0) + ( float(coverage) / abs( float(exon[0])-float(exon[1]) )  )
                break

    for x,y in cov_dict.iteritems():
        record = "%s %.4f\n" % (str(x),float(y))
        c.write(record)
    c.close()
    o.close()
    f.close()

def cov_binner():
    cov_bins = {}
    depths = []
    o = open(covBins, 'w')

    f = open(cov_out)
    depth_max = 0
    for line in f:
        x = line.split(' ')
        depth = float(x[1].strip())
        depths.append(depth)
        if depth > depth_max:
            depth_max = depth

    binS = int(binSize)
    binNum = int(depth_max/binSize)
    for x in range(binNum+1):
        cov_bins[binS] = 0
        binS += binSize
    current_bin = binSize
    for x in np.sort(depths):
        while x > current_bin:
            current_bin += binSize
        cov_bins[current_bin] += 1
    for x,y in cov_bins.iteritems():
        record = "%s %.4f\n" % (x,float(y))
        o.write(record)


#------------------------------------------------------------------------------#

print "Reading in gene locations"
dict_gen()
print "Calculating coverage..."
chr = coverage()
print "Binning Exon Coverages!"
cov_binner()
print "Complete!"
