def recurse(self, positions, history=None):
    from utils import Tree
    left, right, top, bottom = positions
        if history is None:
            history = Tree()
        # right block empty, want to now align this segment
        if top == bottom:
            curr_alignment = [self.v[left:right], '-' * (right - left)]
            return curr_alignment
        # left block empty, want to now align this segment
        elif left == right:
            curr_alignment = ['-' * (bottom - top), self.w[top:bottom]]
            return curr_alignment
        else:
            # get middle node and edge and shift with current source
            node, edge = [
                    tuple(map(sum, zip(_, [left, top])))
                if _ is not None else None
                for _ in middle_edge(
                    self.v[left:right],
                    self.w[top:bottom],
                    self.sigma,
                    self.score_dict)]

            # add middle edge to BST...a BST is necessary for situations where
            # the middle edge is (None, None), i.e. when len(v)==len(w)==1;
            # for these situations, edge is defined by history of traversal
            if (node, edge) == (None, None):
                node, edge = (left, top), (right, bottom)

            # representation of middle node in v and w
            curr_alignment = [
                ['-', self.v[node[0]]][edge[1] - node[1]],
                ['-', self.w[node[1]]][edge[0] - node[0]]
            ]

            # shift edge with current alignment if there was an indel
            if '-' in curr_alignment:
                _node = node
                edge, _node = [tuple(__[_] - 1 if curr_alignment[_] == '-' else __[_] + 1
                                     for _ in range(2)) for __ in (edge, _node)]
                temp = history.insert([_node, edge]).pop(edge)
                history = Tree(elem for elem in history).extend(temp)
            else:
                # add middle edge to tree
                history.extend([node, edge])

            # generate next curr_alignment from left and right subtrees
            l_alignment = recurse(self,
                                  [left, node[0], top, node[1]], history)

            r_alignment = recurse(self,
                                  [edge[0], right, edge[1], bottom], history)

            # create alignment from curr_alignments
            return [l_alignment[i] + curr_alignment[i] + r_alignment[i] for i in range(2)]

    self.v, self.w = recurse(self, [0, len(self.v), 0, len(self.w)])