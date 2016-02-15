# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "anthony"
__date__ = "$Oct 27, 2015 2:55:01 PM$"

import sys
import os.path
from lxml import etree as et
from BedEntry import BedEntry

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

def enter_seqvar(root, acc, seqtype, chromosome, addedInfo, loci, score, geneId, transcriptId, seq):
    entry = et.SubElement(root, UP+'entry', dataset="Ensembl")

    accession = et.SubElement(entry, UP+'accession')
    accession.text = acc

    #Set the names
    name = et.SubElement(entry, UP+'name')
    name.text = seqtype

    fullName = et.SubElement(et.SubElement(et.SubElement(entry, UP+'protein'), UP+'recommendedName'), UP+'fullName')
    fullName.text = addedInfo

    gene = et.SubElement(entry, UP+'gene')
    geneName1 = et.SubElement(gene, UP+'name', type="coords")
    if loci: geneName1.text = chromosome + ":" + loci
    else: geneName1.text = chromosome

    organism = et.SubElement(entry, UP+'organism')
    organismName1 = et.SubElement(organism, UP+'name', type="scientific")
    organismName2 = et.SubElement(organism, UP+'name', type="common")
    organismName1.text = "Homo sapiens"
    organismName2.text = "Human"

    proteinExist = et.SubElement(entry, UP+'proteinExistence', type="evidence at transcript level")
    if score: et.SubElement(proteinExist, UP+'depth', reads=str(score))

    sequence = et.SubElement(entry, UP+'sequence', version="1", fragment="single")
    sequence.text = seq

def generate_tryptic_peps(pep_seq):
    tryptic_peps = []
    k_fragments = pep_seq.split('K')
    for i, k_fragment in enumerate(k_fragments):
        if i != len(k_fragments) - 1: kr_fragments = (k_fragment + 'K').split('R')
        else: kr_fragments = k_fragment.split('R')
        for j, kr_fragment in enumerate(kr_fragments):
            if j != len(kr_fragments) - 1: tryptic_peps.append(kr_fragment + 'R')
            else: tryptic_peps.append(kr_fragment)
    while tryptic_peps and not tryptic_peps[-1]: tryptic_peps.pop() #remove empty elements at end of list
    return tryptic_peps

def check_range(geneModel, chrom, strand, start, end):
    for index in range(start,end + 1):
        if index in geneModel[chrom,strand]: 
            return geneModel[chrom,strand][index]
    
def scan_bed(root, geneModel, splice_bed, nsjDepthCut, refName, scoreName):
    linect = 0
    if splice_bed != None:
        try:
            splice_bed = os.path.abspath(splice_bed)
            linect = sum(1 for line in open(splice_bed))
            splice_bed = open(splice_bed, 'r')
        except Exception, e:
            print >> sys.stderr, "failed: %s" % e
            exit(2)  
        
    try:
        for i, line in enumerate(splice_bed):
            if i % 1000 == 0: print "splice_bed line " + str(i) + " of " + str(linect)
            if line.startswith('track'): continue
            entry = BedEntry(line)
            
            #evaluate for depth cutoff
            if entry.score <= nsjDepthCut: continue

            #Rest of "Scan bed file"
            translations = entry.get_filtered_translations()
            for i,tx in enumerate(translations):
                if tx:
                    (chromStart, chromEnd, translation) = tx

                    #Extended positions
                    ex1Left = entry.chromStart
                    ex2Left = entry.chromStart + entry.blockStarts[1]
                    ex1Right = ex1Left + entry.blockSizes[0] - 1
                    ex2Right = ex2Left + entry.blockSizes[1] - 1
                            
                    trypFragLeft, trypFragRight = chromStart, chromStart - 1
                    for peptide in generate_tryptic_peps(translation):
                        trypFragRight += len(peptide) * 3
                        if trypFragRight > ex1Right and trypFragRight < ex2Left: #exceeds the first exon 
                            trypFragRight = trypFragRight - ex1Right + ex2Left - 1 #set to position in exon2
                        if trypFragLeft < ex1Right and trypFragRight > ex1Right: #Contains the junction, so look up type in gene model, e.g. exon1::novel;exon2:FOX1:protein_coding, else
                            ex1Type = check_range(geneModel, entry.chrom, entry.strand, trypFragLeft, ex1Right)
                            ex2Type = check_range(geneModel, entry.chrom, entry.strand, ex2Left, trypFragRight)
                            if not ex1Type: ex1Type = ('','novel')
                            if not ex2Type: ex2Type = ('','novel')

                            frame_name = '_%s' % (i + 1)
                            acc = entry.name + frame_name
                            chromosome = "chromosome:%s:%s:%s" % (refName, entry.chrom, entry.strand)
                            exons = "exons:%s:%s;%s:%s" % (ex1Type[0], ex1Type[1], ex2Type[0], ex2Type[1])
                            loci = "%s:%s;%s:%s" % (trypFragLeft, ex1Right, ex2Left, trypFragRight)
                            score = entry.score if scoreName == 'depth' else ''
                            enter_seqvar(root, acc, 'pep:splice', chromosome, exons, loci, score, '', '', peptide)

                        trypFragLeft += len(peptide) * 3
                        if trypFragLeft > ex1Right and trypFragLeft < ex2Left:
                            trypFragLeft = trypFragLeft - ex1Right + ex2Left - 1
        splice_bed.close()

    except Exception, e:
        print >> sys.stderr, "failed: Error reading %s - %s" % (splice_bed if splice_bed else 'stdin', e)
    

