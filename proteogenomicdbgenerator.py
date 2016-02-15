# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "anthony"
__date__ = "$Oct 27, 2015 2:41:36 PM$"

import sys
import os.path
import optparse
import refparse
import novelsplices
import variantcalls

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

def write_xml_to_fasta(outFile, root):
    entryct = len(root)
    if outFile != None:
        for i, entry in enumerate(root):
            header = ">"
            if entry.tag == UP+'copyright': continue
            if i % 1000 == 0: print "writing entry " + str(i) + " of " + str(entryct)
            if entry.get('dataset') == 'Ensembl':
                accession, name, geneInfo, chromosome, geneId, transcriptId, addedInfo = None, None, None, None, None, None, None
                accession = entry.find(UP+'accession').text.strip('>')
                name = entry.find(UP + 'name').text
                geneInfo = entry.find(UP+'gene')
                if geneInfo != None: 
                    geneInfo.getiterator(UP+'name')
                    for item in geneInfo:
                        if item.get('type') == 'coords': chromosome = item.text
                        if item.get('type') == 'primary':
                            geneId = 'gene:' + item.text
                            transcriptId = 'transcript:' + item.find(UP+'transcript').text
                addedInfo = entry.find(UP + 'protein').find(UP + 'recommendedName').find(UP + 'fullName').text
                score = entry.find(UP + 'proteinExistence').find(UP + 'depth')
                score = 'depth:' + score.get('reads') if score != None else ''
                headerInfo = [x for x in [accession, name, chromosome, geneId, transcriptId, addedInfo] if x != None] #remove None elements
                header += ' '.join(headerInfo)
            else:
                if entry.get('dataset') == 'Swiss-Prot': database = 'sp'
                else: database = 'tr'
                accession = entry.find(UP + 'accession').text
                name = entry.find(UP + 'name').text
                accession = '|'.join([database, accession, name])
                organism = entry.find(UP + 'organism')
                if organism != None: organism = 'OS=' + organism.find(UP+'name').text
                geneName = entry.find(UP+ 'gene')
                if geneName != None: geneName = 'GN=' + geneName.find(UP+'name').text
                headerInfo = [x for x in [accession, organism, geneName] if x != None] 
                header += ' '.join(headerInfo)
            header = header.replace('\n','').replace('\r','')
            sequence = entry.find(UP+'sequence').text.replace('\n','').replace('\r','')
            outFile.write(header + '\n')
            outFile.write(sequence + '\n')
        outFile.close()
    else: print >> sys.stderr, "write_xml_to_fasta failed: no out-file found"
    
def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
      #I/O
    parser.add_option( '-x', '--reference_xml', dest='reference_xml', help='The reference UniProt-XML file.' )
    parser.add_option( '-p', '--protein_fasta', dest='protein_fasta', help='Reference protein FASTA file.')
    parser.add_option( '-g', '--gene_model', dest='gene_model', default=None, help='GTF gene model file')
    parser.add_option( '-v', '--snpeff_vcf', dest='snpeff_vcf', help='The input snpeff vcf file with HGVS annotations (else read from stdin)' )
    parser.add_option( '-b', '--splice_bed', dest='splice_bed', help='BED file (tophat junctions.bed) with sequence column added' )
    parser.add_option( '-o', '--output', dest='output', help='Output file path. Outputs database in XML format unless --output-fasta is selected.' )
    parser.add_option( '-z', '--output_fasta', dest='output_fasta', action='store_true', default=False, help='Output a FASTA-format database. Place path for file after the --output flag.')
      #Peptide sequence construction
    parser.add_option( '-l', '--leading_aa_num', dest='leading_aa_num', type='int', default=0, help='leading number of AAs to output for SAV peptides' )
    parser.add_option( '-t', '--trailing_aa_num', dest='trailing_aa_num', type='int', default=0, help='trailing number of AAs to output for SAV peptides' )
      #Filtering parameters
    parser.add_option( '-D', '--nsj_depth_cutoff', dest='nsj_depth_cutoff', type='int', default=0, help='Keep only NSJs found with above this depth (DP=# field)' )
    parser.add_option( '-E', '--snv_depth_cutoff', dest='snv_depth_cutoff', type='int', default=0, help='Keep only SNVs found with above this depth (DP=# field)' )
      #Simple entry
    parser.add_option( '-Q', '--bed_score_name', dest='bed_score_name', default="depth", help='include in the NSJ ID line score_name:score, default: "depth."'  )
    parser.add_option( '-R', '--reference', dest='reference', default=None, help='Genome Reference Name for NSJ ID location '  )
    (options, args) = parser.parse_args()
    
    #Protein FASTA
    protein_fasta = refparse.protein_fasta(options.protein_fasta)
    
    #Reference XML/FASTA
    ensembl, uniprot = None, None
    if options.protein_fasta != None: ensembl = refparse.read_fasta_to_xml(options.protein_fasta)
    if options.reference_xml != None: uniprot = refparse.read_and_condense_xml(options.reference_xml)
    
    #Gene model
    geneModel = refparse.gene_model(options.gene_model)

    #Output
    outFile = None
    if options.output == None: outFile = sys.stdout
    else:
        try:
            outFile = os.path.abspath(options.output)
            outFile = open(outFile, 'w')
        except Exception, e:
            print >> sys.stderr, "failed: %s" % e
            exit(3)
    
    root = ensembl.getroot()
    novelsplices.scan_bed(root, geneModel, options.splice_bed, options.nsj_depth_cutoff, options.reference, options.bed_score_name)
    variantcalls.scan_vcf(root, options.snpeff_vcf, protein_fasta, options.snv_depth_cutoff, options.leading_aa_num, options.trailing_aa_num)
    unify_xml(ensembl, uniprot)

    if not options.output_fasta: ensembl.write(outFile, pretty_print=True)
    else: write_xml_to_fasta(outFile, root)

if __name__ == "__main__" : __main__()

