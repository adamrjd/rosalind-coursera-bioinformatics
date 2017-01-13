using System;
using System.Linq;
using System.Collections.Generic;

namespace Coursera
{
    public class LinearSpaceAlignment
    {
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
                        if (v[x - 1] == w[y - 1]) score = s[y - 1, 0] + this.score_dict[this.v[x - 1]][this.w[y - 1]];
                        else score = new int[2] { s[y, 0], s[y - 1, 1] }.Max();

                        s[y, 0] = s[y, 1]; s[y, 1] = score;
                    }
                }
            }
            return s;
        }
        public int[,] MiddleEdge(string v, string w)
        {
            if (w.Length <= 1) return null;
            string c = v;
            v = w;
            w = v;
            int icol = (int)Math.Floor((double)(v.Length + 1) / 2);

            int[,] s = this.LinearScoring(v, w);
            var query = from int[] row
                        in this.LinearScoring(
                            string.Join("", v.Substring(0, icol).Reverse()),
                            string.Join("", w.Reverse()))
                        select row;
            int[][] r_s = query.Reverse().ToArray();

            int[,] edge = new int[2, 2];
            int max_score = -int.MaxValue;

            for (int y = 0; y < icol + 1; y++)
            {
                int node = s[y, 1] + r_s[y][0];
                if (node > max_score)
                {
                    edge[0, 0] = y; edge[0, 1] = icol;
                    max_score = node;
                }
            }

            // if score minimized on diagonal, get diagonal..else get right
            if (r_s[edge[0, 0]][1] < s[edge[0, 0], 1] + r_s[edge[0, 0]][0])
            {
                edge[1, 0] = edge[0, 0] + 1; edge[1, 1] = edge[0, 1] + 1;
            }
            else
            {
                edge[1, 0] = edge[0, 0] + 1; edge[1, 1] = edge[0, 1];
            }

            return edge;
        }

        public void GlobalAlignment()
        {
            string[] global_align(int sigma, Dictionary<char, Dictionary<char, int>> score_dict, string v, string w)
            {
                return new string[2] { "", "" };
            }
            string[] recurse(int top, int bottom, int left, int right)
            {
                if (left == right) return new string[2] { this.v.Substring(top, top + bottom), string.Join("", "-".Take(bottom - top)) };
                else if (top == bottom) return new string[2] { string.Join("", "-".Take(right - left)), this.w.Substring(left, left + right) };
                else if (bottom - top == 1 || right - left == 1) return global_align(this.sigma, this.score_dict, this.v.Substring(top, top + bottom), this.w);
                else
                {
                    int[,] middle_edge = this.MiddleEdge(this.v.Substring(top, top + bottom), this.w.Substring(left, left + right));
                    int[] middle_node = new int[2] { middle_edge[0, 0] + top, middle_edge[0, 1] + left };
                    int[] next_node = new int[2] { middle_edge[1, 0] + top, middle_edge[1, 1] + left };

                    string[] current = new string[2] {
                        new string[2] {"-", this.v.Substring(middle_node[0] % this.v.Length, 1)} [ next_node[0] - middle_node[0] ],
                        new string[2] {"-", this.w.Substring(middle_node[1] % this.w.Length)} [ next_node[1] - middle_node[1] ]};

                    string[] A = recurse(top, middle_node[0], left, middle_node[1]);
                    string[] B = recurse(next_node[0], bottom, next_node[1], right);
                    return new string[2] { A[0] + current[0] + B[0], A[1] + current[1] + B[1] };
                }
            }
            var temp = recurse(0, this.v.Length, 0, this.w.Length);
            this.v_aligned = temp[0]; this.w_aligned = temp[1];
            this.score = this.v_aligned
                .Zip(this.w_aligned,
                    (a, b) => new List<char> { a, b }
                        .Contains('-') ? -this.sigma : this.score_dict[a][b])
                            .Sum();
        }
    }
}