# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

def testvcfscan():
import refparse, variantcalls
p = refparse.protein_fasta('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/Homo_sapiens.GRCh37.73.pep.all.fasta')
d = refparse.read_and_condense_xml('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/SwissProt150105ConciseCondensed.xml')
variantcalls.scan_vcf(d.getroot(), '/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/HepG2SnpEff.vcf', p, 0, 33,33)

def testbedscan():
import refparse, novelsplices
d = refparse.read_and_condense_xml('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/SwissProt150105ConciseCondensed.xml')
g = refparse.gene_model('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/Galaxy10-[Homo_sapiens.GRCh37_canon.73.gtf].gtf')
novelsplices.scan_bed(d.getroot(), g, '/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/HepG2_ExtractGenomicDNA.bed', 2, '', 'depth')

def testensemblxml():
import refparse
p = refparse.protein_fasta('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/Homo_sapiens.GRCh37.73.pep.all.fasta')
d = refparse.read_and_condense_xml('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/SwissProt150105ConciseCondensed.xml')
refparse.ensembl_xml(d.getroot(), p)

def testfastatoxml():
import refparse
d = refparse.read_fasta_to_xml('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/Homo_sapiens.GRCh37.73.pep.all.fasta')

def testwritetofasta():
HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'
import refparse, proteogenomicdbgenerator
d = refparse.read_fasta_to_xml('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/Homo_sapiens.GRCh37.73.pep.all.fasta')
proteogenomicdbgenerator.write_xml_to_fasta(open('../testWriteToFasta.fasta','w'),d.getroot())

def testwritetofasta_fromswissprot():
import refparse, proteogenomicdbgenerator
d = refparse.read_and_condense_xml('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/uniprot-proteome%3AUP000005640.xml')
proteogenomicdbgenerator.write_xml_to_fasta(open('../testWriteToFasta.fasta','w'),d.getroot())

def testcondensexml():
import refparse
d = refparse.read_and_condense_xml('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/uniprot-proteome%3AUP000005640.xml')
d.write('../testing/newcondensedxml.xml',pretty_print=True)

def testunified():
import refparse
p = refparse.read_fasta_to_xml('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/Homo_sapiens.GRCh37.73.pep.all.fasta')
d = refparse.read_and_condense_xml('/home/anthony/NetBeansProjects/ProteogenomicDBGenerator/testing/uniprot-proteome%3AUP000005640.xml')
refparse.unify_xml(p, d)
p.write(open('../testunify.xml','w'),pretty_print=True)

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'
uniprot_root = d.getroot()
uniprot_root.remove(uniprot_root.find(UP+'copyright'))
uniprot_seqs = [x.find(UP+'sequence').text.replace('\n','').replace('\r','') for x in uniprot_root]
len(ensembl_seqs) #84744
len(ensembl_root) #84744
len(set(ensembl_seqs)) #84744
len(uniprot_seqs) #70071
len(uniprot_root) #70071
len(set(uniprot_seqs)) #69948
len(p[1]) #84744
len(set(ensembl_seqs) & set(uniprot_seqs)) #58659
sameLen = 0
for x in uniprot_only:
    for y in ensembl_only:
        if len(x) == len(y): sameLen += 1
sameLen #386907
len(uniprot_only) #11289
len(ensembl_only) #26085



