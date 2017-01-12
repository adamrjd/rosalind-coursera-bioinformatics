using System;
using System.IO;
using System.Linq;
using System.Collections.Generic;

namespace Coursera
{
    public class Program
    {
        public static void Main()
        {
            string path = Directory.GetCurrentDirectory() + "\\text.txt"; //Windows
            //string path = Directory.GetCurrentDirectory() + "//text.txt"; //UNIX
            FileStream stream = new FileStream(path, FileMode.Open);
            StreamReader file = new StreamReader(stream);
            string[] data = file.ReadToEnd().Split(new char[2] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            int n = int.Parse(data[0]);
            List<List<int>> matrix = new Graph(data.ToList().GetRange(1, data.Length - 1)).DistanceBetweenLeaves(n);
            foreach (List<int> row in matrix) System.Console.WriteLine(string.Join(" ", row.Select(x => x.ToString())));
        }
    }
}