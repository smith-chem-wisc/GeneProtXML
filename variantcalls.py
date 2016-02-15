## Amino_Acid_Change notations
# G528R
# p.Gly528Arg/c.1582G>C
#
# This is the current format of the EFF entry:
# EFF=missense(MODERATE|MISSENSE|Ggg/Cgg|G528R|802|SCNN1D|protein_coding|CODING|ENST00000379116|12|1);OICR=(ENST00000379116|1808)
# If this becomes variable, will need to dynamically pattern this on the defintion in the vcf header:
##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank | Genotype_Number [ | ERRORS | WARNINGS ] )' ">

__author__ = "anthony"
__date__ = "$Oct 27, 2015 2:55:09 PM$"

import sys
import os.path
import re
import sys
from lxml import etree as et
import refparse

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'
  
def scan_vcf(root, snpeff_vcf, protein_fasta, snvDepthCut, leading_aas, trailing_aas):
    linect = 0
    if snpeff_vcf != None:
        try:
            snpeff_vcf = os.path.abspath(snpeff_vcf)
            linect = sum(1 for line in open(snpeff_vcf))
            snpeff_vcf = open(snpeff_vcf, 'r')
        except Exception, e:
            print >> sys.stderr, "failed: %s" % e
            exit(2)     
    aa_change_regex = '([A-Z])(\d+)([A-Z])' # G528R
    aa_hgvs_regex = 'p\.([A-Z][a-z][a-z])(\d+)([A-Z][a-z][a-z])(/c\.(\d+)([ACGTN])>([ACGTN]))' # p.Gly528Arg/c.1582G>C
    aa_abbrev_dict = refparse.aa_abbrev_dict()
    vcf_header = [] # Save VCF file header, not currently used
    try:
        for i, line in enumerate(snpeff_vcf):
            if i % 1000 == 0: print "snpeff_vcf line " + str(i) + " of " + str(linect)
            ## print >> sys.stderr, "%d: %s\n" % (linenum,line)
            if line.startswith('##'):
                vcf_header.append(line)
                # May need to check SnpEff version in the header, the EFF info changed between versions 2 and 3
                ##SnpEffVersion
            elif line.startswith('#CHROM'): continue
            else:
                fields = line.split('\t')
                (chrom,pos,id,ref,alts,qual,filter,info) = fields[0:8]
                qual = float(qual)
                for info_item in info.split(';'):
                    try:
                        if info_item.find('=') < 0: continue
                        (key,val) = info_item.split('=',1)
                        if key == 'DP' and int(val) <= snvDepthCut: continue # filter out the entries that have poor coverage AC150310
                        if key == 'EFF':
                            for effect in val.split(','):
                                (eff,effs) = effect.rstrip(')').split('(')                                
                                if eff not in ['NON_SYNONYMOUS_CODING','MISSENSE','missense_variant']: continue #updated for snpeff.4.0 with the inclusion of MISSENSE
                                
                                (impact, functional_class, codon_change, aa_change, aa_len, gene_name, biotype, coding, transcript, exon) = effs.split('|')[0:10]
                                if transcript:
                                    aa_pos = None # 1-based position
                                    alt_aa = '_'
                                    sav = aa_change
                                    m = re.match(aa_change_regex,aa_change) # parse aa_change, and get AA change position and alternate Animo Acid
                                    if m:
                                        aa_pos = int(m.groups()[1])
                                        alt_aa = m.groups()[2]
                                    else:
                                        m = re.match(aa_hgvs_regex, aa_change)
                                        if m:
                                            aa_pos = int(m.groups()[1])
                                            ref_aa = aa_abbrev_dict[m.groups()[0]]
                                            alt_aa = aa_abbrev_dict[m.groups()[2]]
                                            sav = "%s%d%s" % (ref_aa, aa_pos, alt_aa)
                                    if not aa_pos: continue

                                    # get AA sequence
                                    aa_offset, start_pos, end_pos, alt_seq = (aa_pos - 1), 0, 0, ""
                                    
                                    try:
                                        i = 0
                                        for item in protein_fasta[0]:
                                            if item.find(transcript) >= 0: break
                                            i += 1
                                        (pep_id, pep_seq) = (protein_fasta[0][i], protein_fasta[1][i])
                                    except:
                                        print "error finding transcript " + transcript + " in " + line
                                    if not pep_seq: continue
                                    start_pos = max(aa_offset - leading_aas, 0) if leading_aas else 0
                                    end_pos = min(aa_offset + trailing_aas + 1, len(pep_seq)) if trailing_aas else len(pep_seq)
                                    # transform sequence
                                    alt_seq = pep_seq[start_pos:aa_offset] + alt_aa + pep_seq[aa_offset+1:end_pos]

                                    #>ENSP00000317992 pep:sav chromosome:GRCh37:1:879584:894670:-1 gene:ENSG00000188976 transcript:ENST00000327044 gene_biotype:protein_coding transcript_biotype:protein_coding snv_location:1:888659 codon_change:Atc/Gtc sav:I300V
                                    pep_id = re.sub('pep:[a-z]*', 'pep:sap', pep_id)
                                    hdr = ">%s snv_location:%s:%s codon_change:%s sav:%s\n" % (pep_id, chrom, pos, codon_change, sav)
                                    fields = hdr.split(' ')
                                    enter_seqvar(root, fields[0][1:], 'pep:sav', fields[2], ' '.join(fields[5:]), '', '', fields[3][5:], fields[4][11:], alt_seq)
                                    
                    except Exception, e:
                        print >> sys.stderr, "failed: %s" % e
    except Exception, e:
        print >> sys.stderr, "failed: %s" % e
        exit(1)
