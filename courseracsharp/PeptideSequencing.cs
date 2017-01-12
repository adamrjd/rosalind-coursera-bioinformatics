using System.Collections.Generic;
using System.Linq;

namespace Coursera
{
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