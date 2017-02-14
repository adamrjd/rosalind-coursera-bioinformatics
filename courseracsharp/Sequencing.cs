using System.Collections.Generic;
using System.Linq;
using System;

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

    public class PeptideSequencing
    {
        static MassData lookup = new MassData();
        List<List<int>> fragment(List<int> peptide, bool cyclic)
        {
            List<List<int>> fragments = new List<List<int>>();
            if (cyclic)
            {
                for (int i = 0; i < peptide.Count; i++)
                {
                    List<int> pep = peptide.GetRange(i, peptide.Count - i).ToList();
                    // too lazy to extend List<T> for addition
                    foreach (int p in peptide.GetRange(0, i).ToList()) pep.Add(p);
                    for (int j = 1; j < peptide.Count; j++)
                    {
                        fragments.Add(pep.GetRange(0, j));
                    }
                }
                fragments.Add(peptide);
            }
            else
            {
                for (int i = 0; i <= peptide.Count; i++) for (int j = 1; j <= peptide.Count; j++)
                    {
                        fragments.Add(peptide.GetRange(i, peptide.Count - i));
                    }
            }
            return fragments;
        }
        public List<int> theoretical_spectra(List<int> peptide, bool cyclic)
        {
            List<int> theory = new List<int> { 0 };
            foreach (List<int> frag in fragment(peptide, cyclic)) theory.Add(frag.Sum());
            theory.Sort();
            return theory;
        }
        public int spectra_score(List<int> peptide, List<int> spectra, bool cyclic)
        {
            /*Scores peptide sequence against theoretical and experimental spectra by counts
            comparison. If mass of the peptide is too large, returns -1.*/

            //all unique peaks in experimental spectra
            List<int> theoretical = theoretical_spectra(peptide, cyclic);
            HashSet<int> peaks = new HashSet<int>(spectra);
            foreach (int peak in theoretical) peaks.Add(peak);
            //base case
            if (spectra[spectra.Count - 1] < theoretical[theoretical.Count - 1]) return -1;
            //dicts with multiplicities for speed
            Dictionary<int, int> theoretical_mult = new Dictionary<int, int>();
            Dictionary<int, int> spectra_mult = new Dictionary<int, int>();
            foreach (int peak in theoretical)
            {
                if (theoretical_mult.ContainsKey(peak)) theoretical_mult[peak] += 1;
                else theoretical_mult[peak] = 1;
            }
            foreach (int peak in spectra)
            {
                if (spectra_mult.ContainsKey(peak)) spectra_mult[peak] += 1;
                else spectra_mult[peak] = 1;
            }
            //scoring
            int score = 0;
            int[] selection = new int[2];
            foreach (int peak in peaks) score += Enumerable.Min(new int[2] { spectra_mult.ContainsKey(peak) ? spectra_mult[peak] : 0
                , theoretical_mult.ContainsKey(peak) ? theoretical_mult[peak] : 0 });
            return score;
        }
        public string LeaderboardSequencing(List<int> spectra, int N, List<int> aa, bool cyclic = false)
        {
            /*Scores sequences and adds them to a leaderboard, which is trimmed down to size N 
            after every growth cycle. One growth cycle is an appending of every sequence in the
            leaderboard with scoring; growth stops when all possible sequences are too large.
            Returns first candidate from set of highest scoring sequences*/

            Dictionary<int, HashSet<List<int>>> sequences = new Dictionary<int, HashSet<List<int>>>();
            Dictionary<int, HashSet<List<int>>> growths = new Dictionary<int, HashSet<List<int>>>();
            HashSet<List<int>> best_sequences = new HashSet<List<int>>();
            int max_score = 0;
            foreach (int a in aa)
            {
                // seed sequences with initial set of amino acids and scores
                int score = spectra_score(new List<int> { a }, spectra, cyclic);
                if (sequences.ContainsKey(score)) sequences[score].Add(new List<int> { a });
                else sequences[score] = new HashSet<List<int>> { new List<int> { a } };
            }
            while (sequences.Count > 0)
            {
                growths.Clear();
                // grow sequences by single amino acid iff mass of growth is in spectra
                foreach (HashSet<List<int>> vals in sequences.Values) foreach (List<int> seq in vals)
                    {
                        foreach (int a in aa)
                        {
                            List<int> growth = seq.ToList(); // want unique instance, so call function
                            growth.Add(a);
                            int score = spectra_score(growth, spectra, cyclic);
                            if (growths.ContainsKey(score)) growths[score].Add(growth);
                            else growths[score] = new HashSet<List<int>> { growth };
                        }
                    }

                // clear sequences; break if no growths; capture best sequences
                sequences.Clear();
                if (growths.Keys.Max() > max_score)
                {
                    max_score = growths.Keys.Max();
                    // unique instance
                    best_sequences = new HashSet<List<int>>();
                    foreach (List<int> temp in growths[max_score]) best_sequences.Add(temp.ToList());
                }
                int total = growths.Values.Select(x => x.Count).Sum();

                // top N filter with ties
                while (sequences.Values.Select(x => x.Count).Sum() < N && growths.Keys.Count > 0)
                {
                    if (growths.Keys.Count == 0 || growths.Keys.Max() == -1) break;
                    int current_score = growths.Keys.Max();
                    if (growths[current_score].Count == 0)
                    {
                        growths.Remove(current_score);

                        current_score = growths.Keys.Max();
                    }
                    foreach (List<int> spec in growths[current_score])
                    {
                        if (sequences.ContainsKey(current_score)) sequences[current_score].Add(spec);
                        else
                        {
                            sequences[current_score] = new HashSet<List<int>>();
                            sequences[current_score].Add(spec);
                        }
                    }
                    growths.Remove(current_score);
                }
            }
            // first best sequence in string format
            return string.Join("-", best_sequences.First());
        }
        Dictionary<int, int> SpectralConvolutionWithMultiplicities(List<int> spectra)
        {
            Dictionary<int, int> convolution = new Dictionary<int, int>();
            foreach (int i in Enumerable.Range(0, spectra.Count))
            {
                foreach (int mass in spectra.GetRange(0, i))
                {
                    int diff = spectra.ElementAt(i) - mass;
                    if (57 <= diff && diff <= 200)
                    {
                        if (convolution.ContainsKey(diff)) convolution[diff] += 1;
                        else convolution[diff] = 1;
                    }
                }
            }
            return convolution;
        }
        public string ConvolutionCyclopeptideSequencing(int M, int N, List<int> spectra)
        {
            /*Creates a spectral convolution from an experimental spectra, which is then exported
            to the LeaderboardSequencing function to return the first best candidate using scored
            peptide sequences that are grown cyclicly.*/
            Dictionary<int, int> convolution = SpectralConvolutionWithMultiplicities(spectra);
            List<int> convolution_spectra = new List<int>();
            foreach (KeyValuePair<int, int> kvp in convolution.Select(x => x).OrderByDescending(x => x.Value))
            {
                // capture top M multiplicities with ties
                convolution_spectra.Add(kvp.Key);
                if (convolution_spectra.Count == M) break;
            }
            return LeaderboardSequencing(spectra, N, convolution_spectra, true);
        }
    }
}