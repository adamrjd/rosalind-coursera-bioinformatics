using System;
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
            string cwd = Directory.GetCurrentDirectory();
            string path = cwd + "//masstable.txt";
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
    public class Score
    {
        public Dictionary<char, Dictionary<char, int>> NEGATIVE_IDENTITY;
        public Dictionary<char, Dictionary<char, int>> INVERSE_IDENTITY;
        public Dictionary<char, Dictionary<char, int>> IDENTITY;
        public Dictionary<char, Dictionary<char, int>> ALL_ONES;
        public Dictionary<char, Dictionary<char, int>> BLOSUM62;
        public Dictionary<char, Dictionary<char, int>> PAM250;
        public Score()
        {
            string cwd = Directory.GetCurrentDirectory();
            string datdir = "\\datfiles\\";
            string path = cwd + datdir + "//masstable.txt";
            FileStream stream = new FileStream(path, FileMode.Open);
            StreamReader file = new StreamReader(stream);

            // pass streamreader to dictify function and return Dict<char,char> using txt file names as arguments
            Func<string, Dictionary<char, Dictionary<char, int>>> get_dictionary = delegate (string f)
            { return this.dictify(new StreamReader(new FileStream(cwd + datdir + f, FileMode.Open))); };

            this.NEGATIVE_IDENTITY = get_dictionary("NEGATIVE_IDENTITY.txt");
            this.INVERSE_IDENTITY = get_dictionary("INVERSE_IDENTITY.txt");
            this.IDENTITY = get_dictionary("IDENTITY.txt");
            this.ALL_ONES = get_dictionary("ALL_ONES.txt");
            // biologically significant matrices
            this.BLOSUM62 = get_dictionary("BLOSUM62.txt");
            this.PAM250 = get_dictionary("PAM250.txt");
        }
        Dictionary<char, Dictionary<char, int>> dictify(StreamReader datafile)
        {
            Dictionary<char, Dictionary<char, int>> dict = new Dictionary<char, Dictionary<char, int>>();
            string[][] data = datafile.ReadToEnd().Split(new char[] { '\n' }, StringSplitOptions.RemoveEmptyEntries)
                .Select(x => x.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries)).ToArray();
            for (int i = 0; i < data[0].Count(); i++)
            {
                for (int j = 1; j < data.Count(); j++)
                {
                    Dictionary<char, int> temp;
                    if (dict.ContainsKey(data[0][i][0])) temp = dict[data[0][i][0]];
                    else temp = new Dictionary<char, int>();

                    temp.Add(data[j][0][0], int.Parse(data[i + 1][j]));
                    dict[data[0][i][0]] = temp;
                }
            }
            return dict;
        }
    }
}