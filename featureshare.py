# module name: featureshare
# main program: samplespecificdbgenerator

__author__ = "Anthony Cesnik & Michael Knippen"
__date__ = "$Mar 31, 2016 8:00:00 PM$"

from lxml import etree as et
import refparse
import numpy as np

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None : HTML_NS, "xsi" : XSI_NS}
UP = '{'+HTML_NS+'}'

def share_all_features(e_features, e_featureless):
    features = e_features.findall(UP+'feature')
    seq = e_featureless.find(UP+'sequence')
    for feature in features:
        seq.addprevious(feature)

def share_entry_features(entry, ptm_entry, ptm_positions):
    features = ptm_entry.findall(UP+'feature')
    seq = entry.find(UP+'sequence')
    for feature in features:
        if feature.get('type') != 'modified residue':
            seq.addprevious(feature)
        else:
            for location in feature:
                for position in location:
                    if int(position.get('position')) not in (ptm_positions):
                        seq.addprevious(feature)

def aa_num_dict():
    aa_abbrevs = refparse.aa_abbrev_dict().values()
    return {aa_abbrevs[aa_abbrev] : i + 1 for i, aa_abbrev in enumerate(aa_abbrevs)}

def seq_to_num(seq, aa_num_dict):
    seq_as_numbers = np.zeros((1, len(str(seq))))
    for i, ch in enumerate(seq):
        if ch in aa_num_dict():
            seq_as_numbers[0, i] = aa_num_dict[ch]
        else:
            seq_as_numbers[0, i] = 100
    return seq_as_numbers
