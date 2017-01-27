using System;
using System.IO;

namespace Coursera
{
    public class Program
    {
        public static void Main()
        {
            string cwd = Directory.GetCurrentDirectory();
            FileStream stream = new FileStream(cwd + "\\input.txt", FileMode.Open);
            StreamReader file = new StreamReader(stream);
            var data = file.ReadToEnd().Split(new char[] { '\n', '\r' }, System.StringSplitOptions.RemoveEmptyEntries);
            string str1 = data[0];
            string str2 = data[1];

            var answer = new Coursera.LinearSpaceAlignment(
                str1,
                str2,
                5,
                new Coursera.Score().BLOSUM62);
            answer.GlobalAlignment();
            Console.WriteLine(answer.score);
            Console.WriteLine(answer.v_aligned);
            Console.WriteLine(answer.w_aligned);

        }
    }
}