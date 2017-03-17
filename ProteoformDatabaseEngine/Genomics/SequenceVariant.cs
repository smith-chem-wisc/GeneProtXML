using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Genomics
{
    public class SequenceVariant : ChromosomeSegment
    {
        public double qual { get; set; }
        public string reference { get; set; }
        public string alternate { get; set; }
        public double allele_frequency { get; set; }
        public int depth { get; set; }
        public new int length { get { return Math.Max(reference.Length, alternate.Length); } }
        public SequenceVariant(Chromosome chrom, int position, string id, string reference, string alternate, double qual, double allele_frequency, int depth) 
            : base(id, chrom, "+", position, position + Math.Max(reference.Length, alternate.Length) - 1, null, null)
        {
            this.qual = qual;
            this.reference = reference;
            this.alternate = alternate;
            this.allele_frequency = allele_frequency;
            this.depth = depth;
        }

        public override string ToString()
        {
            return String.Join(":", new string[] { this.chrom.name, this.start.ToString(), this.end.ToString(), this.reference, this.alternate, this.allele_frequency.ToString() });
        }
    }

    public class SNV : SequenceVariant
    {
        public SNV(Chromosome chrom, int position, string id, string reference, string alternate, double qual, double allele_frequency, int depth)
            : base(chrom, position, id, reference, alternate, qual, allele_frequency, depth)
        { }

        public bool is_missense()
        {
            return false;
        }
        public bool is_start_gain()
        {
            return false;
        }
        public bool is_stop_loss()
        {
            return false;
        }

        public bool parse_aa_change(string aa_change)
        {
        //aa_abbrev_dict = amino_acids_3to1
        //aa_change_regex = '([A-Z])(\d+)([A-Z])'  # G528R
        //aa_hgvs_regex = 'p\.([A-Z][a-z][a-z])(\d+)([A-Z][a-z][a-z])(/c\.(\d+)([ACGTN])>([ACGTN]))'  # p.Gly528Arg/c.1582G>C
        //aa_pos = None  # 1-based position
        //ref_aa, alt_aa = '_', '_'
        //m = re.match(aa_change_regex, aa_change)  # parse aa_change, and get AA change position and alternate Animo Acid
        //if m:
        //    aa_pos = int(m.groups()[1])
        //    ref_aa = m.groups()[0]
        //    alt_aa = m.groups()[2]
        //else:
        //    m = re.match(aa_hgvs_regex, aa_change)
        //    if m:
        //        aa_pos = int(m.groups()[1])
        //        ref_aa = aa_abbrev_dict[m.groups()[0]]
        //        alt_aa = aa_abbrev_dict[m.groups()[2]]
        //return aa_pos, ref_aa, alt_aa
            return false;
        }
    }

    public class Indel : SequenceVariant
    {
        public Indel(Chromosome chrom, int position, string id, string reference, string alternate, double qual, double allele_frequency, int depth)
            : base(chrom, position, id, reference, alternate, qual, allele_frequency, depth)
        { }
    }
    public class Insertion : Indel
    {
        public Insertion(Chromosome chrom, int position, string id, string reference, string alternate, double qual, double allele_frequency, int depth)
            : base(chrom, position, id, reference, alternate, qual, allele_frequency, depth)
        { }
    }
    public class Deletion : Indel
    {
        public Deletion(Chromosome chrom, int position, string id, string reference, string alternate, double qual, double allele_frequency, int depth)
            : base(chrom, position, id, reference, alternate, qual, allele_frequency, depth)
        { }
    }
}
