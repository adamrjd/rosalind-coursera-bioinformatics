using System.Collections.Generic;

namespace Extensions
{
    public static class Generic
    {
        public static string[] Split(string s, string chars)
        {
            // sloppy extension b/c i hate \r\n
            List<string> array = new List<string>();
            int i = 0;
            int start = 0;
            int len = chars.Length;
            while (i <= s.Length - len)
            {
                if (s.Substring(i, len) == chars)
                {
                    array.Add(s.Substring(start, i - start));
                    i += len;
                    start = i;
                }
                else i++;
            }
            if (start < s.Length - len) array.Add(s.Substring(start, s.Length - start));
            return array.ToArray();
        }
    }
}