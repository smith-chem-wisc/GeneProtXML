# module name: indel
# main program: samplespecificdbgenerator

__author__ = "Michael Knippen"
__date__ = "$Mar 30, 2016 2:55:03 PM$"

import collections as col
import numpy as np
from lxml import etree as et

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

def append_indel(entry, count, result, seq_len, fasta_len):

    # Find the sequence tag to add the modified residue before
    tmp = entry.find(UP+"sequence")
    index = entry[entry.index(tmp) - 1]

    #amino acid and location of the indel
    aa = result[0]
    pos = result[1]

    #case where an amino acid has been deleted
    if (seq_len > fasta_len):
        feature = et.Element(UP+"feature", type="sequence variant", description="Found in FASTA", id='"'+str(count)+'"')
        var = et.Element(UP+"deletion")
        var.text = aa
        feature.append(var)

    #case where an amino acid has been added
    else:
        feature = et.Element(UP+"feature", type="sequence variant", description="Found in FASTA", id='"'+str(count)+'"')
        var = et.Element(UP+"insertion")
        var.text = aa
        feature.append(var)

    location = et.Element(UP + "location")
    position = et.Element(UP + "position", position='"' + pos + '"')
    location.append(position)

    feature.append(location)

    index.addnext(feature)

def append_seqvar(entry, count, tmp):
    #Find the sequence tag to add the modified residue before
    tmp = entry.find(UP+"sequence")
    index = entry[entry.index(tmp)-1]

    count += 1
    x = np.nonzero(tmp)[1]
    original_seq = "'" + seq[x] + "'" #TODO: seq is an unresolved reference
    variation_seq = "'" + str(fasta_dict[key].seq)[x] + "'" #TODO: declare fasta_dict or bring it in as a praramter

    #Create a lxml element parts
    feature = et.Element(UP+"feature", type="sequence variant", description="Found in FASTA", id='"'+str(count)+'"')
    original = et.Element(UP+"original")
    original.text = original_seq
    variant = et.Element(UP+"variation")
    variant.text = variation_seq
    location = et.Element(UP+"location")
    position = et.Element(UP+"position", position='"'+str(x)+'"')
    location.append(position)

    #Create the entry
    feature.append(original)
    feature.append(variant)
    feature.append(location)

    #Add seqvar to the db
    index.addnext(feature)

#Process Sequence
def count_aa(seq):
    d = col.Counter()
    for char in seq:
        if char in d.keys(): d[char] += 1
        else: d[char] = 1
    return d

#TODO: make some boolean variables to clarify what you're evaluating in these if/elses
def swap_seq(seq1, seq2):    
    swap = False #should the values be switched, due to negative and 0 counter values removed
    if len(seq1 - seq2) == 1:
        swap = False
    elif len(seq2 - seq1) == 1:
        swap = True
    else:
        swap = None
    return swap

#TODO: make variable names more specific
#Find Larger Sequence    
def indel(seq1, seq2):
    c1 = count_aa(seq1)
    c2 = count_aa(seq2)
    swap = swap_seq(c1,c2)
    if swap == True: #sub
        c3 = c2-c1  
    elif swap == False: 
        c3 = c1-c2
    else: c3 = None
    return c3
    
def indel_location(seq1, seq2, aa):
    larger = seq2
    smaller = seq1
    if len(seq1)>len(seq2):
       larger = seq1
       smaller = seq2
    run = True
    length = len(larger)
    start = 0
    while run == True:        
        pos = larger.find(aa, start, length)
        if larger[pos] == smaller[pos]:
            start = pos+1
        else:
            run == False
            return (aa,pos)



#TODO: make this method name more descriptive
def main(str1,str2): #0 based indexing used
    result = indel(str1, str2)
    if result == None: return None
    else: return indel_location(str1, str2, result.keys()[0])
