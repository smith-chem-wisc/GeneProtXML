using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using Proteomics;
using GenomicsData;
using System.IO;
using Proteogenomics;

namespace Test
{
    [TestFixture]
    public class TestSequenceSimilarityMethods
    {
        [Test]
        public void test_modification_transfer_exact_sequence_match()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            Dictionary<string, Modification> un;
            List<Protein> ok = ProteinDbLoader.LoadProteinDb(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml2.xml"), false, nice, false, out un);
            List<Protein> destination = new List<Protein> { new Protein("MKTCYYELLGVETHASDLELKKAYRKKALQYHPDKNPDNVEEATQKFAVIRAAYEVLSDPQERAWYDSHKEQILNDTPPSTDDYYDYEVDATVTGVTTDELLLFFNSALYTKIDNSAAGIYQIAGKIFAKLAKDEILSGKRLGKFSEYQDDVFEQDINSIGYLKACDNFINKTDKLLYPLFGYSPTDYEYLKHFYKTWSAFNTLKSFSWKDEYMYSKNYDRRTKREVNRRNEKARQQARNEYNKTVKRFVVFIKKLDKRMKEGAKIAEEQRKLKEQQRKNELNNRRKFGNDNNDEEKFHLQSWQTVKEENWDELEKVYDNFGEFENSKNDKEGEVLIYECFICNKTFKSEKQLKNHINTKLHKKNMEEIRKEMEEENITLGLDNLSDLEKFDSADESVKEKEDIDLQALQAELAEIERKLAESSSEDESEDDNLNIEMDIEVEDVSSDENVHVNTKNKKKRKKKKKAKVDTETEESESFDDTKDKRSNELDDLLASLGDKGLQTDDDEDWSTKAKKKKGKQPKKNSKSTKSTPSLSTLPSSMSPTSAIEVCTTCGESFDSRNKLFNHVKIAGHAAVKNVVKRKKVKTKRI",
                "", new Dictionary<int, List<Modification>>(), null, null, null, "", "", false, false, new List<GoTerm>()) };

            Assert.AreEqual(ok[0].BaseSequence, destination[0].BaseSequence);
            List<Protein> new_proteins = SequenceSimilarity.TransferModifications(ok, destination);

            Assert.AreEqual(ok[0].OneBasedBeginPositions, new_proteins[0].OneBasedBeginPositions);
            Assert.AreEqual(ok[0].OneBasedEndPositions, new_proteins[0].OneBasedEndPositions);
            Assert.AreEqual(ok[0].OneBasedPossibleLocalizedModifications, new_proteins[0].OneBasedPossibleLocalizedModifications);
            Assert.True(new_proteins[0].OneBasedPossibleLocalizedModifications.Keys.Count == 2);
            Assert.AreEqual(ok[0].GoTerms.Count, new_proteins[0].GoTerms.Count);
            Assert.AreEqual(ok[0].BaseSequence, new_proteins[0].BaseSequence);
            Assert.AreEqual(1, new_proteins.Count);
        }
    }
}
