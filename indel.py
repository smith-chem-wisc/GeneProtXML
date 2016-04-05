# module name: indel
# main program: samplespecificdbgenerator

__author__ = "Michael Knippen"
__date__ = "$Mar 30, 2016 2:55:03 PM$"

import collections as col
import numpy as np
import featureshare
from lxml import etree as et

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

def append_indel(entry, uniprot, result):
    x = []
    for item in result:
        x.append(int(item[1]))
    featureshare.share_entry_features(entry,uniprot, x)

def append_seqvar(entry, uniprot, comparison):
    x = np.nonzero(comparison)
    x = x[1]
    x = np.squeeze(np.asarray(x))
    featureshare.share_entry_features(entry, uniprot, x)

#Process Sequence
def count_aa(seq):
    d = col.Counter()
    for char in seq:
        if char in d.keys(): d[char] += 1
        else: d[char] = 1
    return d


#TODO: make variable names more specific DONE:
#Find Larger Sequence    
def indel(ens, uni):
    ens_ct = count_aa(ens)
    uni_ct = count_aa(uni)
    ens_del = ens_ct-uni_ct
    if ens_del.__len__() == 0: ens_del = None
    ens_in = uni_ct-ens_ct
    if ens_in.__len__() == 0: ens_in = None
    indels = (ens_del, ens_in)
    return indels
    
def indel_location(seq1, seq2, aa, search_position, insertion):
    larger = seq2
    smaller = seq1
    if len(seq1)>len(seq2):
       larger = seq1
       smaller = seq2
    run = True
    length = len(larger)
    if search_position == 0: start = 0
    else: start = search_position
    while run == True:
        if insertion:
            pos = larger.find(aa, start, length)
            if larger[pos] == smaller[pos]: start = pos + 1
            else:
                run == False
                return (aa, pos)
        else:
            pos = smaller.find(aa, start, length)
            if larger[pos] == smaller[pos]: start = pos+1
            else:
               run == False
               return (aa,pos)

def update_results(indel_locations, deletions):
    tmp_indels = sorted(indel_locations, key=lambda x:x[1])
    print tmp_indels
    updated_indels = []
    in_del = col.Counter({'in':0,'del':0})
    for i in range(tmp_indels.__len__()):
        in_del_counter = in_del['in']-in_del['del']
        update_position = list(tmp_indels[i])
        update_position[1] = str(int(update_position[1])+in_del_counter)
        updated_indels.append(tuple(update_position))
    print updated_indels
    return updated_indels

def append_indel_info(ens_seq, uni_seq, result, indel_loc):
    num_indels = 0
    for key in result[0].keys(): num_indels += result[0][key]
    for key in result[1].keys(): num_indels -= result[1][key]
    indel_locations = indel_loc
    no_repeated_aa_indel=True
    insertion = True
    for item in result:
        for key in item.keys():
            if item[key] > 1: no_repeated_aa_indel = False
    if no_repeated_aa_indel:
        for key in result[0].keys(): indel_locations.append(indel_location(ens_seq,uni_seq,key,0, insertion))
        insertion = False
        for key in result[1].keys(): indel_locations.append(indel_location(ens_seq,uni_seq,key,0,insertion))
        indel_locations = update_results(indel_locations, result[1])
        return indel_locations
    else: return None #TODO: Allow for repeated indels

def reinsert_deletion(sequence, indel_locations):
    seq = list(sequence)
    for i in range(len(indel_locations)): seq.insert(int(indel_locations[i][1]),indel_locations[i][0])
    seq = ''.join(seq)
    return seq

def delete_insertion(sequence, indel_locations):
    seq = list(sequence)
    for i in range(len(indel_locations)): seq.__delitem__(int(indel_locations[i][1]))
    seq = ''.join(seq)
    return seq

def validate_seqs(ens, uni, indel_locations, result):
    insertions = []
    deletions = []
    for item in indel_locations:
        if item[0] in result[0].keys(): deletions.append(item)
        elif item[0] in result[1].keys(): insertions.append(item)
    uni_seq = uni
    uni_seq = delete_insertion(uni_seq,insertions)
    uni_seq = reinsert_deletion(uni_seq,deletions)
    if ens == uni_seq: return True
    else: return False

#TODO: make this method name more descriptive
def get_indel_locations(ens,uni): #0 based indexing used
    result = indel(ens, uni)
    if result == None: return None
    else:
        indel_locations = []
        indel_locations = append_indel_info(ens, uni, result, indel_locations)
        if validate_seqs(ens,uni,indel_locations,result): return indel_locations
        else: return None
