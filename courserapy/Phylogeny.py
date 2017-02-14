'''
Module for Coursera Bioinformatics IV
Created by: Adam Dornford
Written in: Py35
'''


class Graph(object):
    '''Creates graphs in python collections from genetic data'''
    from Utils import Utils
    adjlist = Utils.nested_dict()
    leaves = list()
    matrix = None
    phylogeny = None

    def __init__(self, data, graph_type="weighted_graph"):
        if graph_type == "weighted_graph":
            try:
                # data already consumed
                assert isinstance(data, dict)
                self.adjlist = data
            except KeyError:
                # data = list("0->1:2","1->2:3",..."row_n->etc->etc")
                from re import split
                for edge in data:
                    _ = list(filter(lambda _: _ != '',
                                    split(r'[:\->]', edge)))
                    _ = _ if len(_) == 3 else _.append(
                        1)  # if no weight, add 1
                    self.adjlist[_[0]][_[1]] = _[2]
                vals = [_ for __ in [val.keys() for val in [self.adjlist[k] for k in self.adjlist]]
                        for _ in __]
                self.leaves = [leaf for leaf in set(
                    vals) if vals.count(leaf) <= 1]
        elif graph_type == "distance_matrix":
            # data = list("0 1 2 3","1 0 1 2",...,"row_n")
            self.matrix = [[int(x) for x in row.split(" ")] for row in data]
            for i, row in enumerate(data):
                for j, weight in enumerate(row.split(" ")):
                    self.adjlist[int(i)][int(j)] = int(weight)
            for _ in range(len(self.adjlist)):
                self.leaves.append(_)

    def adjlist_dict_to_strlist(self):
        if self.phylogeny is not None:
            self.adjlist = self.phylogeny
        if len(self.adjlist.keys()) != 0:
            adj_str_ls = list()
            while len(self.adjlist.keys()) > 0:
                item = self.adjlist.popitem()
                k = item[0]
                for v in item[1]:
                    adj_str_ls.append(
                        str(k) + '->' + str(v) + ':' + str(item[1][v]))
            return adj_str_ls
        return

    def distance_between_leaves(self, n):
        def Dijkstra(source, sink):
            nodes = list(self.adjlist)
            dist_dict = {v: float('inf') for v in nodes}
            dist_dict[source] = 0
            prev_node = dict()
            while len(nodes) > 0:
                u = [node for node in nodes if dist_dict[node] == min(
                    [dist_dict[v] for v in dist_dict if v in nodes])][0]
                nodes.remove(u)
                for v in self.adjlist[u]:
                    alt = dist_dict[u] + self.adjlist[u][v]
                    if alt < dist_dict[v]:
                        dist_dict[v] = alt
                        prev_node[v] = u
            return dist_dict[sink]
        if self.matrix != None:
            self.matrix = [[0 for _ in range(n)] for __ in range(n)]
            for source in self.leaves:
                for sink in self.leaves:
                    if source > sink:
                        self.matrix[source][sink] = Dijkstra(source, sink)
                        self.matrix[sink][source] = self.matrix[source][sink]
        return self.matrix

    def limb_length(self, n, j):
        d = lambda x, y: self.matrix[x][y]
        limb_list = list()
        for i in range(n):
            if i != j:
                for k in range(n):
                    if i != k != j:
                        limb_list.append(((d(i, j) + d(j, k) - d(i, k)) / 2))
        return int(min(limb_list))

    def additive_phylogeny(self, n):
        from copy import deepcopy

        # lambdas
        dist_theorem = lambda i, j, k: (
            self.adjlist[i][j] + self.adjlist[j][k] - self.adjlist[i][k]) / 2
        limb_length = lambda node: int(min([dist_theorem(i, node - 1, k)
                                            for i in range(len(self.leaves))
                                            for k in range(len(self.leaves))
                                            if i != node - 1 != k and i != k]))

        if n == 2:
            # base case
            return self

        graph = Graph(deepcopy(self.adjlist))
        length = limb_length(n)
        #___ bald tree ___
        for _ in graph.adjlist.keys():
            if _ != n - 1:
                graph.adjlist[n - 1][_] -= length
                graph.adjlist[_][n - 1] -= length

        i, k = None, None
        for _, __ in zip(list(sorted(graph.adjlist.keys())),
                         list(sorted(graph.adjlist.keys()))[::-1]):
            if graph.adjlist[_][__] == graph.adjlist[_][n - 1] + graph.adjlist[n - 1][__]:
                i, k = _, __
                break

        # new limb distance to find
        x = graph.adjlist[_][n - 1]

        #___ trim tree ___
        for _ in graph.adjlist.keys():
            graph.adjlist[_].pop(n - 1)
        graph.adjlist.pop(n - 1)
        graph.leaves.remove(n - 1)

        # recurse to get new adjlist
        T = graph.additive_phylogeny(n - 1)

        #___ add limb to T ___

        v = max(T.leaves) + 1

        # create new limb
        __ = 0
        for _ in range(v):
            __ += T.adjlist[0][_]
            if x - __ == 0:
                # attach node as leaf
                (T.adjlist[_][v],
                 T.adjlist[v][_],
                 T.adjlist[v][v]) = length, length, 0
                for node in list(i for i in range(v - 1) if i != _):
                    dist = int(dist_theorem(v, node, _))
                    T.adjlist[v][node] = dist
                self.leaves.append(max(self.leaves) + 1)
            elif x - __ < 0:
                # insert node as internal node

        print('yay')
        # return new adjlist
        self.phylogeny = T
        return T


if __name__ == "__main__":
    matrix_print = lambda m: [
        print('\t'.join([str(_) for _ in row])) for row in m]
    with open('input.txt', 'r') as f:
        N, *INPUT_DATA = [l.strip() for l in f.readlines()]
    N = int(N)
    M = Graph(INPUT_DATA, "distance_matrix")
    M.additive_phylogeny(N)
    for _ in M.adjlist_dict_to_strlist():
        print(_)
