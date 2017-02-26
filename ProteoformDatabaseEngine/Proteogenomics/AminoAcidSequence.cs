using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ProteoformDatabaseEngine
{
    public class AminoAcidSequence
    {
        public string name { get; set; }
        public string original_sequence { get; set; }
        public string sequence { get; set; }
        public Dictionary<int, Modification> modifications { get; set; } = new Dictionary<int, Modification>();
        public List<SequenceVariant> sequence_variants { get; set; } = new List<SequenceVariant>();
        public List<Modification> other_features { get; set; } = new List<Modification>();
        public Transcript transcript { get; set; }
        public int trans_start { get; set; }
        public int trans_end { get; set; }

        public int length { get { return sequence.Length; } }
        public override string ToString() { return this.sequence; }

        public AminoAcidSequence(string name, string sequence, Transcript transcript, int trans_start, int trans_end)
        {
            this.name = name;
            this.original_sequence = sequence;
            this.sequence = sequence;
            this.transcript = transcript;
            this.trans_start = trans_start;
            this.trans_end = trans_end;
        }
        public AminoAcidSequence(string name, string sequence) : this(name, sequence, null, -1, -1)
        { }
            
        public List<AminoAcidSequence> tryptic_digestion()
        {
            List<AminoAcidSequence> return_seqs = new List<AminoAcidSequence>();
            List<string> tryptic_peps = new List<string>();
            string[] k_fragments = this.sequence.Split(new char[] { 'K' });
            for (int i = 0; i < k_fragments.Length; i++)
            {
                string k_frag = k_fragments[i];
                string[] kr_fragments;
                if (i != k_fragments.Length - 1)
                    kr_fragments = (k_frag + "K").Split(new char[] { 'R' });
                else 
                    kr_fragments = k_frag.Split(new char[] { 'R' });
                for (int j = 0; j < kr_fragments.Length; j++)
                {
                    string kr_frag = kr_fragments[j];
                    if (j != kr_fragments.Length - 1)
                        tryptic_peps.Add(kr_frag + "R");
                    else
                        tryptic_peps.Add(kr_frag);
                }
            }
            tryptic_peps = tryptic_peps.Where(x => !string.IsNullOrEmpty(x)).ToList();

            int tryp_start = this.trans_start;
            int tryp_end = this.trans_end;
            for (int i = 0; i < tryptic_peps.Count; i++)
            {
                string tryptic_pep = tryptic_peps[i];
                tryp_end += tryptic_pep.Length * 3;
                return_seqs.Add(new AminoAcidSequence(this.name + "_" + i.ToString(), tryptic_pep, this.transcript, tryp_start, tryp_end);
                tryp_start += tryptic_pep.Length * 3;
            }
            return return_seqs;
        }

        public void trim(int start, int end)
        {
            //todo: check that it can actually be trimmed
            this.sequence = this.sequence.Substring(start, end - start);
            this.trans_start += start * 3;
            this.trans_end -= (this.length - end) * 3;
        }

        public void methionine_cleavage()
        {
            if (this.sequence[0] == 'M') this.trim(1, this.length);
        }
    }
}
