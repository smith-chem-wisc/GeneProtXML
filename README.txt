Author: Anthony Cesnik
Creation: October 27, 2015
Last updated: October 28, 2015

ProteogenomicDBGenerator is a program that takes RNA sequencing data analysis results and translates them into protein database entries that are appended to a UniProt-XML file. This XML is condensed, meaning that much of the extraneous author information in UniProt-XMLs is removed for faster search times in Morpheus.

ProteogenomicDBGenerator is based on the following papers:
1. Sheynkman, et al. "Discovery and Mass Spectrometric Analysis of Novel Splice-Junction Peptides Using RNA-Seq." Mol Cell Proteomics 2013, 12, 2341-2353.
2. Sheynkman, et al. "Large-scale mass spectrometric detection of variant peptides resulting from nonsynonymous nucleotide differences." J Proteome Research 2014, 13, 228-240.
3. Sheynkman, et al. "Using Galaxy-P to leverage RNA-Seq for the discovery of novel protein variations." BMC Genomics 2014, 15, 9.
4. Cesnik, et al. "Human Proteomic Variation Revealed by Combining RNA-Seq Proteogenomics and Global Post-Translational Modification (G-PTM) Search Strategy." In review.

v0.0.1
October 27, 2015
This program breaks several scripts into modules, translate_bed_sequences.py based on (1) and snpefftopeptides.py based on (2), which were reported in (3). In (4), these were appended to the UniProt-XML to garner modifed and sequence variant peptide identifications. In that work, we also found that filtering using the read depth proved effective, where methods coded into the previous scripts were not. Only the read depth filters were left in this script. Finally, only splice junction peptides containing no stop codon in the first exon are included.

v0.0.2
October 28, 2015
Implemented using a FASTA database as the main reference input and as the output format. Bug fix with read_and_condense_xml where the organism field was being deleted.

