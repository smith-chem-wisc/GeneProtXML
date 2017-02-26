using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ProteoformDatabaseEngine
{
   public  class SpliceJunction
    {
        int intron_start { get; set; }
        int intron_end { get; set; }

        public override bool Equals(object obj)
        {
            return base.Equals(obj);
        }
    }
}
