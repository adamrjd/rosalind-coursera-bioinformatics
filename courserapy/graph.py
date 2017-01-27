'''
Module for Coursera Bioinformatics IV
Created by: Adam Dornford
Written in: Py35
'''


class Graph(object):
    '''Creates graphs in python collections from genetic data'''
    from utils import nested_dict
    adjlist = nested_dict()
    leaves = list()
    matrix = None
    phylogeny = None

    def __init__(self, data, graph_type="weighted_graph"):
        if graph_type == "weighted_graph":
            # data = list("0->1:2","1->2:3",..."row_n->etc->etc")
            from re import split
            for edge in data:
                _ = list(filter(lambda _: _ != '',
                                split(r'[:\->]', edge)))
                _ = _ if len(_) == 3 else _.append(1)  # if no weight, add 1
                self.adjlist[_[0]][_[1]] = _[2]
            vals = [_ for __ in [val.keys() for val in self.adjlist.values()]
                    for _ in __]
            self.leaves = [leaf for leaf in set(vals) if vals.count(leaf) <= 1]
        elif graph_type == "distance_matrix":
            # data = list("0 1 2 3","1 0 1 2",...,"row_n")
            self.matrix = [[int(x) for x in row.split(" ")] for row in data]
            for i, row in enumerate(data):
                for j, weight in enumerate(row.split(" ")):
                    self.adjlist[int(i)][int(j)] = int(weight)
            for _ in range(len(self.adjlist.keys())):
                self.leaves.append(_)

    def adjlist_dict_to_strlist(self):
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
            nodes = list(self.adjlist.keys())
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

    def additive_phylogeny(self, curr_adjlist=None, edges=None, length=None):
        from copy import deepcopy
        dist_theorem = lambda i, j, k: (
            self.adjlist[i][j] + self.adjlist[j][k] - self.adjlist[i][k]) / 2
        limb_length = lambda node: int(min([dist_theorem(i, node, k)
                                            for i in range(len(self.leaves))
                                            for k in range(len(self.leaves))
                                            if i != node != k]))

        if curr_adjlist is None:
            curr_adjlist = deepcopy(self.adjlist)

        # add base case...
        if len(curr_adjlist.keys()) > 2:
            node = min(curr_adjlist.keys())
            length = limb_length(node)
            # bald tree
            for _ in curr_adjlist.keys():
                curr_adjlist[node][_] -= length
                if _ != node:
                    curr_adjlist[_][node] -= length
            # trim tree
            edges = curr_adjlist.pop(node)
            for _ in curr_adjlist.keys():
                curr_adjlist[_].pop(node)
            # recurse
            self.additive_phylogeny(curr_adjlist, edges, length)

        # add things back to tree
        if self.phylogeny is not None:
            # recover results from last recursion
            curr_adjlist = self.phylogeny

        for node in curr_adjlist:

        self.phylogeny = deepcopy(curr_adjlist)
        return

if __name__ == "__main__":
    matrix_print = lambda m: [
        print('\t'.join([str(_) for _ in row])) for row in m]
    with open('input.txt', 'r') as f:
        N, *INPUT_DATA = [l.strip() for l in f.readlines()]
    N = int(N)
    M = Graph(INPUT_DATA, "distance_matrix")
    M.additive_phylogeny()
    for _ in M.adjlist_dict_to_strlist():
        print(_)
