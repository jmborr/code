import os, string, sys
from networkx import *

def bfs(G, source):
    D = {}
    P = [source]
    D[source] = 0
    seen = {}
    queue = [source]
    seen[source] = True
    while queue:
        v = queue.pop(0)
        for w in G.neighbors(v):
            if w not in seen:
                seen[w] = True
                queue.append(w)
                D[w] = (D[v] + G.get_edge(v, w))
                P.append((v, w))
                return D[w], P
    return 0, P


if __name__ == '__main__':
    G = XGraph()
    G.add_edge(1, 2, 0.1)
    G.add_edge(1, 4, 0.3)
    # G.add_edge(1, 3)
    G.add_edge(2, 3, 0.2)
    G.add_edge(2, 4, 0.2)
    G.add_edge(4, 3, 0.1)
    l = bfs(G, 1)
    print l