'''
Module for Coursera Bioinformatics IV
Created by: Adam Dornford
Written in: Python35
'''


class Graph(object):
    '''Creates graphs in python collections from genetic data'''
    from utils import nested_dict
    adjlist = nested_dict()
    leaves = list()
    matrix = None  # sometimes consuming matrix not graph

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
                    if i != j:
                        self.adjlist[int(i)] = {int(j): int(weight)}
        elif graph_type == "genome":
            # data = list("+1","-2",...,"+4")
            pass

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
            dist_dict = {v: float('inf') for v in self.adjlist.keys()}
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

    def additive_phylogeny(self, n, matrix):
        dist_theorem = lambda m, i, j, k: (m[i][j] + m[j][k] - m[i][k]) / 2
        limb_length = lambda n, m: int(min([dist_theorem(m, i, n - 1, k)
                                            for i in range(n)
                                            for k in range(n)
                                            if i != k and k != n - 1 and i != n - 1]))
        node = n + 1
        adj_lst = None
        for i in range(n - 2):
            j = i + 1
            k = j + 1
            curr_adj_lst = dict()
            curr_adj_lst[i] = {node: dist_theorem(matrix, i, j, k)}
            curr_adj_lst[j] = {node: dist_theorem(matrix, j, k, i)}
            curr_adj_lst[k] = {node: dist_theorem(matrix, k, i, j)}
            curr_adj_lst[node] = {_: curr_adj_lst[_] for _ in (i, j, k)}
            if adj_lst != None:
                p_node = node - 1
                edges = adj_lst.pop(p_node)
                # want to call limb length for new limb...then insert into old
                # adjlist
            else:
                adj_lst = dict()
                for item in curr_adj_lst:
                    adj_lst.update(item)
            node += 1
        self.adjlist = adj_lst
        return self.adjlist
if __name__ == "__main__":
    matrix_print = lambda m: [
        print('\t'.join([str(_) for _ in row])) for row in m]
    with open('text.txt', 'r') as f:
        N, *INPUT_DATA = [l.strip() for l in f.readlines()]
    N = int(N)
    M = Graph(INPUT_DATA, "distance_matrix")
    matrix_print(M.additive_phylogeny(N, M.matrix))
