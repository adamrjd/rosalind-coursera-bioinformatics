using System;
using System.Linq;
using System.Collections.Generic;

namespace Coursera
{
    public class LinearSpaceAlignment
    {
        /*Translation from python3 b/c py was too slow for 500+ peptide sequences,
        but then I realized both have a bug somewhere ):*/
        string v;
        string w;
        public string v_aligned;
        public string w_aligned;
        public int score;
        int sigma;
        Dictionary<char, Dictionary<char, int>> score_dict;
        public LinearSpaceAlignment(string v, string w, int sigma, Dictionary<char, Dictionary<char, int>> score_dict)
        {
            this.v = v;
            this.w = w;
            this.sigma = sigma;
            this.score_dict = score_dict;

        }
        int[,] LinearScoring(string v, string w)
        {
            int[,] s = new int[w.Length + 1, 2];
            for (int x = 0; x < v.Length + 1; x++)
            {
                for (int y = 0; y < w.Length + 1; y++)
                {
                    if (y == 0 || x == 0) { s[y, 0] = -this.sigma * y; s[y, 1] = -this.sigma * x; }

                    else
                    {

                        int score;
                        if (v[x - 1] == w[y - 1]) score = s[y - 1, 0] + this.score_dict[v[x - 1]][w[y - 1]];
                        else score = new int[2] { s[y, 0], s[y - 1, 1] }.Max();

                        s[y, 0] = s[y, 1];
                        s[y, 1] = score;
                    }
                }
            }
            return s;
        }
        public int[,] MiddleEdge(string v, string w)
        {
            if (w.Length <= 1) return null;
            var c = v;
            v = w;
            w = c;
            int icol = (int)Math.Floor((double)(v.Length) / 2);

            int[,] s = this.LinearScoring(v.Substring(0, icol), w);
            int[,] r_s = this.LinearScoring(
                            string.Join("", v.Substring(icol, v.Length - icol).Reverse()),
                            string.Join("", w.Reverse()));

            int[,] edge = new int[2, 2];
            int max_score = -int.MaxValue;

            for (int y = 0; y < icol; y++)
            {
                int node = s[y, 1] + r_s[w.Length - 1 - y, 0];
                if (node > max_score)
                {
                    edge[0, 0] = y; edge[0, 1] = icol;
                    max_score = node;
                }
            }

            // if score minimized on diagonal, get diagonal..else get right
            if (r_s[w.Length - 1 - edge[0, 0], 1] < s[edge[0, 0], 1] + r_s[w.Length - 1 - edge[0, 0], 0])
            {
                edge[1, 0] = edge[0, 0] + 1; edge[1, 1] = edge[0, 1] + 1;
            }
            else
            {
                edge[1, 0] = edge[0, 0] + 1; edge[1, 1] = edge[0, 1];
            }

            return edge;
        }
        public string[] nonlinear_global_align(string v, string w)
        {
            int i = v.Length;
            int j = w.Length;
            int[,] s = new int[i + 1, j + 1];
            int[,] backtrack = new int[i + 1, j + 1];

            for (int x = 1; x < i + 1; x++)
            {
                for (int y = 1; y < j + 1; y++)
                {
                    int[] l = new int[3] {
                            s[x - 1, y] - this.sigma,
                            s[x, y - 1] - this.sigma,
                            s[x - 1, y - 1] + this.score_dict[v[x - 1]][w[y - 1]]
                        };
                    s[x, y] = l.Max();
                    backtrack[x, y] = l.ToList().FindIndex(0, 1, elem => elem == s[x, y]);
                }
            }

            this.score = s[i, j];
            while (i * j != 0)
            {
                if (backtrack[i, j] == 0)
                {
                    i -= 1;
                    w = w.Insert(j, "-");
                }
                else if (backtrack[i, j] == 1)
                {
                    j -= 1;
                    v = v.Insert(i, "-");
                }
                else
                {
                    i -= 1;
                    j -= 1;
                }
            }
            for (int _ = 0; _ < i; _++) w = w.Insert(0, "-");
            for (int _ = 0; _ < j; _++) v = v.Insert(0, "-");
            return new string[2] { v, w };
        }
        string[] recurse(int top, int bottom, int left, int right)
        {
            if (left == right) return new string[2] {
                this.v.Substring(top, bottom - top),
                "".PadRight(right - left, '-')
            };
            else if (top == bottom) return new string[2] {
                "".PadRight(right - left, '-'),
                this.w.Substring(left, right - left)
            };
            else if (bottom - top == 1 || right - left == 1) return this.nonlinear_global_align(
                this.v.Substring(top, bottom - top),
                this.w.Substring(left, right - left)
            );
            else
            {
                int[,] middle_edge = this.MiddleEdge(this.v.Substring(top, bottom - top), this.w.Substring(left, right - left));
                int[] middle_node = new int[2] { middle_edge[0, 0] + top, middle_edge[0, 1] + left };
                int[] next_node = new int[2] { middle_edge[1, 0] + top, middle_edge[1, 1] + left };

                string[] current = new string[2] {
                        new string[2] {"-", this.v.Substring(middle_node[0] % this.v.Length, 1)} [ next_node[0] - middle_node[0] ],
                        new string[2] {"-", this.w.Substring(middle_node[1] % this.w.Length, 1)} [ next_node[1] - middle_node[1] ]
                        };

                string[] A = recurse(top, middle_node[0], left, middle_node[1]);
                string[] B = recurse(next_node[0], bottom, next_node[1], right);
                return new string[2] { A[0] + current[0] + B[0], A[1] + current[1] + B[1] };
            }
        }
        public void GlobalAlignment()
        {
            var temp = this.recurse(0, this.v.Length, 0, this.w.Length);

            this.v_aligned = temp[0];
            this.w_aligned = temp[1];

            this.score = this.v_aligned
                .Zip(this.w_aligned,
                    (a, b) => new List<char> { a, b }
                        .Contains('-') ? -this.sigma : this.score_dict[a][b])
                            .Sum();
        }
    }
}