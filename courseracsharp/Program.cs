using System;

namespace Coursera
{
    public class Program
    {
        public static void Main()
        {
            string str2 = "PLEASANTLY";
            string str1 = "MEANLY";
            var answer = new Coursera.LinearSpaceAlignment(str1, str2, 5, new Coursera.Score().BLOSUM62);
            answer.GlobalAlignment();
            Console.WriteLine(answer.score);
            Console.WriteLine(answer.v_aligned);
            Console.WriteLine(answer.w_aligned);
        }
    }
}