using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeneProtXml
{
    class IndividualHuman
    {
        //IUPAC letters
        static string amino_acids = "ACDEFGHIKLMNPQRSTVWY";
        static Dictionary<char, string> amino_acids_1to3 = new Dictionary<char, string>() {
            { 'A', "Ala" }, { 'C', "Cys" }, { 'D', "Asp" }, { 'E', "Glu" }, { 'F', "Phe" }, { 'G', "Gly" }, { 'H', "His" }, { 'I', "Ile" }, { 'K', "Lys" },
            { 'L', "Leu" }, { 'M', "Met" }, { 'N', "Asn" }, { 'P', "Pro" }, { 'Q', "Gln" }, { 'R', "Arg" }, { 'S', "Ser" }, { 'T', "Thr" }, { 'V', "Val" },
            { 'W', "Trp" }, { 'Y', "Tyr" }
        };

        static Dictionary<string, char> amino_acids_3to1 = new Dictionary<string, char>() {
            { "Ala", 'A' }, { "Cys", 'C' }, { "Asp", 'D' }, { "Glu", 'E' }, { "Phe", 'F' }, { "Gly", 'G' }, { "His", 'H' }, { "Ile", 'I' }, { "Lys", 'K' },
            { "Leu", 'L' }, { "Met", 'M' }, { "Asn", 'N' }, { "Pro", 'P' }, { "Gln", 'Q' }, { "Arg", 'R' }, { "Ser", 'S' }, { "Thr", 'T' }, { "Val", 'V' },
            { "Trp", 'W' }, { "Tyr", 'Y' }
        };

        static string ambiguous_dna_letters = "GATCRYWSMKHBVDN";
        static string unambiguous_dna_letters = "GATC";
        static Dictionary<string, string> ambiguous_dna_values = new Dictionary<string, string>() {
            { "A", "A" }, { "C", "C" }, { "G", "G" }, { "T", "T" },
            { "M", "AC" }, { "R", "AG" }, { "W", "AT" }, { "S", "CG" }, { "Y", "CT" }, { "K", "GT" },
            { "V", "ACG" }, { "H", "ACT" }, { "D", "AGT" }, { "B", "CGT" }, { "X", "GATC" }, { "N", "GATC" }
        };
        static Dictionary<char, char> ambiguous_dna_complement = new Dictionary<char, char>() {
            { 'A', 'T' },
            { 'C', 'G' },
            { 'G', 'C' },
            { 'T', 'A' },
            { 'M', 'K' },
            { 'R', 'Y' },
            { 'W', 'W' },
            { 'S', 'S' },
            { 'Y', 'R' },
            { 'K', 'M' },
            { 'V', 'B' },
            { 'H', 'D' },
            { 'D', 'H' },
            { 'B', 'V' },
            { 'X', 'X' },
            { 'N', 'N' }
        };
        static Dictionary<char, char> unambiguous_dna_complement = new Dictionary<char, char>() {
            { 'A', 'T' },
            { 'C', 'G' },
            { 'G', 'C' },
            { 'T', 'A' },
        };
        static List<string> standard_start_codons = new List<string>() { "ATG" };
        static List<string> extended_start_codons = new List<string>() { "TTG", "CTG", "ATG" };
        static List<string> standard_stop_codons = new List<string>() { "TAA", "TAG", "TGA" };
        static List<string> mitochon_start_codons = new List<string>() { "ATT", "ATC", "ATA", "ATG", "GTG" };
        static List<string> mitochon_stop_codons = new List<string>() { "TAA", "TAG", "AGA", "AGG" };
        static char[] standard_base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG".ToCharArray();
        static char[] standard_base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG".ToCharArray();
        static char[] standard_base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG".ToCharArray();
        static char[] standard_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG".ToCharArray();
        static char[] mitochon_acids = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG".ToCharArray();

        static Dictionary<string, char> standard_code = Enumerable.Range(0, standard_acids.Length).ToDictionary(
            i => standard_base1[i].ToString() + standard_base2[i].ToString() + standard_base3[i].ToString(),
            i => standard_acids[i]);
        static Dictionary<string, char> mitochon_code = Enumerable.Range(0, mitochon_acids.Length).ToDictionary(
            i => standard_base1[i].ToString() + standard_base2[i].ToString() + standard_base3[i].ToString(),
            i => mitochon_acids[i]);

    }
}
