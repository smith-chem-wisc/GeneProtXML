using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace UnitTestProject1
{
    [TestClass]
    public class IUPACTests
    {
        [TestMethod]
        public void ListsContainOnlyCanonicalMonomers()
        {
            foreach char c in GeneProtXml.Individual
                Assert.IsTrue
        }
    }
}
