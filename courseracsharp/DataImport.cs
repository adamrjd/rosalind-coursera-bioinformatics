using System.IO;
using System.Collections.Generic;
using System.Linq;

namespace Coursera
{
    public class MassData
    {
        static Dictionary<char, int> acid_to_mass = new Dictionary<char, int>();
        static Dictionary<int, char> mass_to_acid = new Dictionary<int, char>();
        public MassData()
        {
            //string path = Directory.GetCurrentDirectory() + "\\masstable.txt";
            string path = Directory.GetCurrentDirectory() + "//masstable.txt";
            FileStream stream = new FileStream(path, FileMode.Open);
            StreamReader file = new StreamReader(stream);
            parse(file);
        }
        static void parse(StreamReader file)
        {
            foreach (string line in file.ReadToEnd().Split('\n'))
            {
                string[] temp = line.Split(' ');
                acid_to_mass[temp[0][0]] = int.Parse(temp[1]);
            }
            foreach (KeyValuePair<char, int> x in acid_to_mass)
            {
                if (!mass_to_acid.ContainsKey(x.Value)) mass_to_acid[x.Value] = x.Key;
            }
        }
        public int get_mass(char acid) { return acid_to_mass[acid]; }
        public char get_acid(int mass) { return mass_to_acid[mass]; }
        public List<int> masses() { return mass_to_acid.Keys.ToList(); }
        public List<char> acids() { return acid_to_mass.Keys.ToList(); }
    }
}