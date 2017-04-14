class Edge():
    def __init__(self, flow, capacity, u, v):
        self.flow = flow
        self.capacity = capacity
        self.u = u
        self.v = v

class Vertex(): 
    def __init__(self, h, edge_flow):
        self.h = h
        self.edge_flow = edge_flow

class Graph():
    def __init__(self, num_vertices, E):
        self.V = list()
        self.E = E
        for i in range(num_vertices):
            self.V.append(Vertex(0, 0))

    def push(u):
        for e in self.E:
            if e.u == u:
                if e.flow == e.capacity:
                    continue
                if self.V[u] > self.V[e.v].h:
                    flow = min(e.capacity - e.flow, V[u].edge_flow)
                    V[u].edge_flow -= flow
                    V[e.v].edge_flow += flow
                    e.flow += flow
                    j = 0
                    for reverse in self.E:
                        if reverse.v == e.u and reverse.u == e.v:
                            reverse.flow -= flow
                            break
                        j += 1
                    if j == len(self.E):
                        self.E.append(Edge(0, flow, e.v, e.u))
                    return True
        return False

    def relabel(u):
        min_height = float("inf")
        for e in self.E:
            if e.u == u:
                if e.flow == e.capacity:
                    continue
                if self.V[e.v].h < min_height:
                    min_height = self.V[e.v].h
                    self.V[u].h = min_height + 1

g = Graph(5, [(0, 4, 20), (1, 2, 3), (3, 1, 4)])
