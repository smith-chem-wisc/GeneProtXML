using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeneProtXml
{
    public class Chromosome
    {
        public string name { get; set; }
        public string sequence { get; set; }
        public List<Gene> genes { get; set; }
        public List<AminoAcidSequence> amino_acid_sequences { get; set; }
        public int length { get { return sequence.Length; } }
        public override string ToString()
        {
            return this.sequence;
        }

        public bool contains(string gene_name)
        {
            return (from gene in genes select gene.name).Contains(gene_name);
        }

        public void sort()
        {
            //TODO
        }

        public Gene get_gene_by_position(int position)
        {
            return new Gene();
        }
    }
}
