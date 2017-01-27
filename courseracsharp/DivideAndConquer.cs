using System;
using System.Linq;
using System.Collections.Generic;

namespace Coursera
{
    public class LinearSpaceAlignment
    {
        /*Translation from cpy35 b/c py is too slow for large sequences!!!*/
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
        Tuple<Dictionary<Tuple<int, int>, int>, int[]> LinearScoring(string v, string w)
        {
            /*score sequences and only preserve a certain number of rows in scoring matrix*/

            // instantiate and initialize matrices
            Dictionary<Tuple<int, int>, int> s = Enumerable.Range(0, v.Length + 1)
                                        .ToDictionary(x => new Tuple<int, int>(x, 0), x => -this.sigma * x);
            int[] backtrack = Enumerable.Range(0, v.Length + 1)
                                        .Select(x => 1)
                                        .ToArray();

            // fill matrices
            for (int y = 1; y < w.Length + 1; y++)
            {
                for (int x = 0; x < v.Length + 1; x++)
                {
                    if (x == 0) s.Add(new Tuple<int, int>(x, y), -this.sigma * y);

                    else
                    {
                        List<int> scores = new List<int> {
                            s[new Tuple<int, int>(x - 1, y - 1)]
                            + this.score_dict[
                                v.Substring(x - 1, 1)[0]][
                                    w.Substring(y - 1, 1)[0]], //diagonal
                            s[new Tuple<int, int>(x, y - 1)] - this.sigma, //right
                            s[new Tuple<int, int>(x - 1, y)] - this.sigma //down
                        };
                        int score = scores.Max();
                        s.Add(new Tuple<int, int>(x, y), score);
                        backtrack[x] = scores.IndexOf(score);
                    }

                    if (y - 2 >= 0 && v.Length > 2) s.Remove(new Tuple<int, int>(x, y - 2));
                }
            }
            return new Tuple<Dictionary<Tuple<int, int>, int>, int[]>(s, backtrack);
        }
        List<int[]> getEdge(string v, string w)
        {
            int jcol = (int)Math.Floor((double)(w.Length) / 2);

            Dictionary<Tuple<int, int>, int> s = this.LinearScoring(
                            v,
                            w.Substring(0, jcol)).Item1;
            var reversed_s = this.LinearScoring(
                            string.Join("", v.Reverse()),
                            string.Join("", w.Substring(jcol, w.Length - jcol).Reverse()));
            Dictionary<Tuple<int, int>, int> r_s = reversed_s.Item1;
            int[] backtrack = reversed_s.Item2;

            int[] node = new int[2];
            int max_score = -int.MaxValue;

            for (int x = 0; x < (int)Math.Floor((decimal)v.Length / 2) + 1; x++)
            {
                int score = s[new Tuple<int, int>(x, jcol)] + r_s[new Tuple<int, int>(v.Length - x, jcol)];
                if (score > max_score)
                {
                    node[0] = x;
                    node[1] = jcol;
                    max_score = score;
                }
            }

            int dir = backtrack[v.Length - node[0]];
            int[] edge;

            if (dir == 0) edge = new int[2] { node[0] + 1, node[1] + 1 };
            else if (dir == 1) edge = new int[2] { node[0] + 1, node[1] };
            else edge = new int[2] { node[0], node[1] + 1 };

            return new List<int[]> { node, edge };
        }
        public List<int[]> MiddleEdge(string v, string w)
        {
            if (v.Length == 1 && w.Length == 1)
                return new List<int[]> { new int[2] { 0, 0 }, new int[2] { 1, 1 } };
            else if (v.Length <= 2 && w.Length > 1)
                return getEdge(
                    w, v)
                    .Select(p => p.Reverse().ToArray())
                    .ToList();
            else
                return getEdge(v, w);
        }
        string[] recurse(int left, int right, int top, int bottom)
        {
            if (top == bottom)
                return new string[2] {
                this.v.Substring(left, right - left),
                "".PadRight(right - left, '-')
            };
            else if (left == right)
                return new string[2] {
                "".PadRight(bottom - top, '-'),
                this.w.Substring(top, bottom - top)
            };

            else
            {
                List<int[]> middle_edge = this.MiddleEdge(
                    this.v.Substring(left, right - left),
                    this.w.Substring(top, bottom - top));
                int[] node = new int[2] {
                    middle_edge[0][0] + left,
                    middle_edge[0][1] + top };
                int[] edge = new int[2] {
                    middle_edge[1][0] + left,
                    middle_edge[1][1] + top };

                string[] subalignment = new string[2] {
                        new string[2] {"-", this.v.Substring(node[0], 1)} [ edge[1] - node[1] ],
                        new string[2] {"-", this.w.Substring(node[1], 1)} [ edge[0] - node[0] ]
                        };

                // lefthand side of divide and conquer recursive traversal
                string[] A = recurse(left, node[0], top, node[1]);

                // shift edge with current subalignment if there is an indel
                if (subalignment.Contains("-"))
                {
                    for (int i = 0; i < 2; i++)
                    {
                        if (subalignment[i].Contains("-")) edge[i] -= 1;
                        else edge[i] += 1;
                    }
                }

                // righthand side of divide and conquer recursive traversal
                string[] B = recurse(edge[0], right, edge[1], bottom);

                // returns results at node in relative recursive frame
                return new string[2] { A[0] + subalignment[0] + B[0], A[1] + subalignment[1] + B[1] };
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