using System;
using System.Collections.Generic;
using System.Linq;

namespace Coursera
{
    public class Graph
    {
        public Dictionary<int, Dictionary<int, int>> AdjacencyList = new Dictionary<int, Dictionary<int, int>>();
        public HashSet<int> Leaves = new HashSet<int>();
        public Graph(List<string> data)
        {
            this.Consume(data);
            List<int> vals = this.AdjacencyList.Values.SelectMany(x => x.Keys).ToList();
            foreach (int leaf in vals) if (vals.Count(x => x == leaf) <= 1) this.Leaves.Add(leaf);
        }
        void Consume(List<string> data)
        {
            foreach (string edge in data)
            {
                List<int> temp = edge.Split("->:".ToArray(), StringSplitOptions.RemoveEmptyEntries).Select(x => int.Parse(x)).ToList();
                if (!edge.Contains(":")) temp.Add(1);
                if (!this.AdjacencyList.ContainsKey(temp[0]))
                {
                    Dictionary<int, int> kv = new Dictionary<int, int>();
                    kv.Add(temp[1], temp[2]);
                    this.AdjacencyList[temp[0]] = kv;
                }
                else this.AdjacencyList[temp[0]].Add(temp[1], temp[2]);
            }
        }
        int Dijkstra(int source, int sink)
        {
            List<int> Q = this.AdjacencyList.Keys.ToList();
            List<int> nodes = this.AdjacencyList.Keys.ToList();
            Dictionary<int, int> DistanceDict = new Dictionary<int, int>();
            foreach (int node in nodes) if (node != source) DistanceDict.Add(node, int.MaxValue); else DistanceDict.Add(node, 0);
            Dictionary<int, int> PrevNode = new Dictionary<int, int>();

            while (Q.Count > 0)
            {
                int u = Q.Select(x => x)
                    .Where(x => DistanceDict[x] ==
                        Q.Select(y => DistanceDict[y]).Min())
                            .ElementAt(0);
                Q.Remove(u);
                foreach (KeyValuePair<int, int> val in this.AdjacencyList[u])
                {
                    int v = val.Key;
                    int alt = DistanceDict[u] + this.AdjacencyList[u][v];
                    if (alt < DistanceDict[v])
                    {
                        DistanceDict[v] = alt;
                        if (!PrevNode.ContainsKey(v)) PrevNode.Add(v, u);
                        else PrevNode[v] = u;
                    }
                }
            }
            return DistanceDict[sink];

        }
        public List<List<int>> DistanceBetweenLeaves(int n)
        {
            List<List<int>> Matrix = new List<List<int>>();
            foreach (int i in Enumerable.Range(0, n)) Matrix.Add(Enumerable.Repeat(0, n).ToList());
            foreach (int source in this.Leaves)
            {
                foreach (int sink in this.Leaves)
                {
                    if (source > sink)
                    {
                        Matrix[source][sink] = this.Dijkstra(source, sink);
                        Matrix[sink][source] = Matrix[source][sink];
                    }
                }
            }
            return Matrix;
        }
    }
}