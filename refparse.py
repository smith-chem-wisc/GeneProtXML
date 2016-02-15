# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "anthony"
__date__ = "$Oct 27, 2015 2:54:44 PM$"

import sys
import os.path
from lxml import etree as et
import variantcalls

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None : HTML_NS, "xsi" : XSI_NS}
UP = '{'+HTML_NS+'}'

def aa_abbrev_dict():
## dictionary for Amino Acid Abbreviations
    aa_abbrev_dict = dict()
    aa_abbrev_dict['Phe'] = 'F'
    aa_abbrev_dict['Leu'] = 'L'
    aa_abbrev_dict['Ser'] = 'S'
    aa_abbrev_dict['Tyr'] = 'Y'
    aa_abbrev_dict['Cys'] = 'C'
    aa_abbrev_dict['Trp'] = 'W'
    aa_abbrev_dict['Pro'] = 'P'
    aa_abbrev_dict['His'] = 'H'
    aa_abbrev_dict['Gln'] = 'Q'
    aa_abbrev_dict['Arg'] = 'R'
    aa_abbrev_dict['Ile'] = 'I'
    aa_abbrev_dict['Met'] = 'M'
    aa_abbrev_dict['Thr'] = 'T'
    aa_abbrev_dict['Asn'] = 'N'
    aa_abbrev_dict['Lys'] = 'K'
    aa_abbrev_dict['Val'] = 'V'
    aa_abbrev_dict['Ala'] = 'A'
    aa_abbrev_dict['Asp'] = 'D'
    aa_abbrev_dict['Glu'] = 'E'
    aa_abbrev_dict['Gly'] = 'G'
    return aa_abbrev_dict

def read_and_condense_xml(xml_file_name):
    if xml_file_name != None:
        try:
            reference_xml = os.path.abspath(xml_file_name)
            refXml = open(reference_xml, 'r')
            p = et.XMLParser(remove_blank_text=True) #required for pretty additions
            db = et.parse(refXml, p)
            root = db.getroot()
            for entry in root:
                for element in entry:
                    if element.tag not in [UP+'protein',UP+'accession',UP+'name',UP+'gene',UP+'organism',UP+'proteinExistence',UP+'depth',UP+'sequence',UP+'feature',UP+'dbReference']:
                        entry.remove(element)
                    elif element.get('type') != 'modified residue' and element.tag == UP+'feature': entry.remove(element)
                    elif element.get('type') != 'Ensembl' and element.tag == UP+'dbReference': entry.remove(element)
                    elif element.tag == UP+'organism':
                        for field in element:
                            if field.tag != UP+'name':
                                element.remove(field)
                    elif element.tag == UP+'protein':
                        for name in element:
                            if name.tag != UP+'recommendedName':
                                element.remove(name)
                    else:
                      continue
            return db
        except Exception, e:
            print >> sys.stderr, "failed: %s" % e
            exit(2)
    
def add_unified(root, newId, rootIndex, acc, seqtype, chromosome, biotypes, geneId, transcriptId, seq):
    if newId: 
        entry = et.SubElement(root, UP+'entry', dataset="Ensembl")
        accession = et.SubElement(entry, UP+'accession')
        accession.text = str(newId)
    else: entry = root[rootIndex]
    
    accession = et.Element(UP+'accession')
    accession.text = acc
    a = entry.findall(UP+'accession')
    if len(a) > 0: a[-1].addnext(et.Element(UP+'accession'))
    else: entry.append(accession)
    
    name = et.Element(UP+'name')
    name.text = seqtype
    n = entry.findall(UP+'name')
    if len(n) > 0: n[-1].addnext(et.Element(UP+'name'))
    else: entry.append(name)

    if newId:
        fullName = et.SubElement(et.SubElement(et.SubElement(entry, UP+'protein'), UP+'recommendedName'), UP+'fullName')
        fullName.text = 'referenceProtein' + str(newId)
        gene = et.SubElement(entry, UP+'gene')
    else:
        gene = et.Element(UP+'gene')
        entry.findall(UP+'gene')[-1].addnext(gene)
    geneName1 = et.SubElement(gene, UP+'name', type="coords")
    geneName1.text = chromosome
    geneName2 = et.SubElement(gene, UP+'name', type="primary")
    geneName2.text = geneId
    gene_biotype = et.SubElement(gene, UP+'gene_biotype')
    gene_biotype.text = biotypes.split(' ')[0]
    transcript = et.SubElement(gene, UP+'transcript', type="primary")
    transcript.text = transcriptId
    transcript_biotype = et.SubElement(gene, UP+'transcript_biotype')
    transcript_biotype.text = biotypes.split(' ')[1]

    if newId:
        organism = et.SubElement(entry, UP+'organism')
        organismName1 = et.SubElement(organism, UP+'name', type="scientific")
        organismName2 = et.SubElement(organism, UP+'name', type="common")
        organismName1.text = "Homo sapiens"
        organismName2.text = "Human"
        proteinExist = et.SubElement(entry, UP+'proteinExistence', type="evidence at transcript level")
        sequence = et.SubElement(entry, UP+'sequence', version="1", fragment="single")
        sequence.text = seq
    
def ensembl_entry(uniqSeqs, root, acc, seqtype, chromosome, biotypes, geneId, transcriptId, seq):
    found = False
    for i, s in enumerate(uniqSeqs):
        if seq == s:
            found = True
            add_unified(root, '', i, acc, seqtype, chromosome, biotypes, geneId, transcriptId, seq)
            break
    if not found:
        uniqSeqs.append(seq)
        add_unified(root, str(len(uniqSeqs)), -1, acc, seqtype, chromosome, biotypes, geneId, transcriptId, seq)
    
def read_fasta_to_xml(protein_fasta):
    if protein_fasta != None:
        try:
            protein_fasta = os.path.abspath(protein_fasta)
            refFasta = open(protein_fasta, 'r')
            root = et.Element(UP+'uniprot', nsmap=NAMESPACE_MAP)
            db = et.ElementTree(root)
            uniqSeqs = []
            
            seq = ""
            line = refFasta.readline().strip()
            while line != "":
                if line.startswith(">"):
                    line = line.split(' ')
                    acc, seqtype, chromosome = line[0][1:], line[1], line[2]
                    geneId, transcriptId, biotypes = line[3].split(':')[1], line[4].split(':')[1], ' '.join(line[5:])
                    line = refFasta.readline().strip()
                    while not line.startswith(">"):
                        if line == "": break
                        seq += line
                        line = refFasta.readline().strip()
                    if seq.find('*') < 0: ensembl_entry(uniqSeqs, root, acc, seqtype, chromosome, biotypes, geneId, transcriptId, seq.replace('\n','').replace('\r',''))
                    seq = ""
            return db
        except Exception, e:
            print >> sys.stderr, "failed: %s" %e
            exit(2)
    
def unify_xml(ensembl, uniprot):
    #TODO: this might be more efficient if I got all of the nodes up front and then just got the parents...
    ensembl_root, uniprot_root = ensembl.getroot(), uniprot.getroot()
    uniprot_root.remove(uniprot_root.find(UP+'copyright'))
    ensembl_seqs = [x.find(UP+'sequence').text.replace('\n','').replace('\r','') for x in ensembl_root] #lists all sequences without linebreaks
    uniprot_seqs = [x.find(UP+'sequence').text.replace('\n','').replace('\r','') for x in uniprot_root] #lists all sequences without linebreaks

    common_seqs = list(set(ensembl_seqs) & set(uniprot_seqs))
    uniprot_only = list(set(uniprot_seqs) - set(common_seqs))
    ensembl_only = list(set(ensembl_seqs) - set(common_seqs))
    
    for s in common_seqs:
        e_entry, u_entry = None, None
        for i, t in enumerate(ensembl_seqs):
            if s == t:
                e_entry = ensembl_root[i]
                break
        for i, t in enumerate(uniprot_seqs):
            if s == t:
                u_entry = uniprot_root[i]
                break
        mods = u_entry.findall(UP+'feature')
        seq = e_entry.find(UP+'sequence')
        for mod in mods:
            seq.addprevious(mod)
    return ensembl
            
def gene_model(gene_model_file_name):
    if gene_model_file_name != None:
        try:
            geneModelPath = os.path.abspath(gene_model_file_name)
            geneModelFile = open(geneModelPath, 'r')
            geneModel = {} # list of [start,end] entries to each [chrom,strand]
            for line in geneModelFile:
                line = line.split('\t')
                chromStart,chromEnd,chrom,strand = int(line[3]),int(line[4]),line[0],line[6]
                attribsTmp = line[8].split('; ')
                attribs = {} #This is what's consistent across Ensembl gene models...
                for item in attribsTmp:
                    item = item.split(' ')
                    #print item
                    attribs[item[0]] = item[1][1:-1]
                #print attribsTmp
                if 'transcript_biotype' not in attribs: attribs['transcript_biotype'] = line[1] #JJ's canonical gtf
                try:
                    for i in range(chromStart,chromEnd+1): geneModel[chrom,strand][i] = [attribs['transcript_name'],attribs['transcript_biotype']]
                except KeyError:
                    geneModel[chrom,strand] = {}
                    for i in range(chromStart,chromEnd+1): geneModel[chrom,strand][i] = [attribs['transcript_name'],attribs['transcript_biotype']]
            return geneModel
        except Exception, e:
            print >> sys.stderr, "failed: %s" % e
            exit(2)
         