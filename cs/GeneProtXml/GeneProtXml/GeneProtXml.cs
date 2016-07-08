using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using ManyConsole;

namespace GeneProtXml
{
    class GeneProtXml
    {
        //[STAThread]
        static void Main(string[] args)
        {
            var commands = GetCommands();
            if (args.Length > 0)
            {
                ConsoleCommandDispatcher.DispatchCommand(commands, args, Console.Out);
            }
            else
            {
                //Application.EnableVisualStyles();
                //Application.SetCompatibleTextRenderingDefault(false);
                //Application.Run(new Form1());
                ConsoleCommandDispatcher.DispatchCommand(commands, new string[] { "user_interface" }, Console.Out);
            }
        }

        public static IEnumerable<ConsoleCommand> GetCommands()
        {
            return ConsoleCommandDispatcher.FindCommandsInSameAssemblyAs(typeof(GeneProtXml));
        }
    }
}
