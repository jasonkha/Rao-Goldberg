import numpy as np
import math
import sys
from collections import deque

class Edge():
    def __init__(self, v, flow, capacity, reverse):
        self.v = v
        self.flow = flow
        self.capacity = capacity
        self.reverse = reverse

class Graph():
    def __init__(self, num_vertices, E):
        self.num_vertices = num_vertices
        self.adj = dict()
        self.level = [0] * num_vertices
        self.num_edges = 0
        self.max_capacity = 0
        self.edges = dict()
        self.resid_cap = self.get_capacities()

    def add_edge(self, u, v, cap):
        try:
            u_queue = self.adj[u]
        except KeyError:
            self.adj[u] = deque([])
            u_queue = self.adj[u]

        try:
            v_queue = self.adj[v]
        except KeyError:
            self.adj[v] = deque([])
            v_queue = self.adj[v]

        forward = Edge(v, 0, cap, len(self.adj[v]))
        reverse = Edge(u, 0, 0, len(self.adj[u]))
        self.edges[(u, v)] = forward

        u_vertex_numbers = [u.v for u in u_queue]
        v_vertex_numbers = [v.v for v in v_queue]
        if v not in u_vertex_numbers:
            u_queue.append(forward)
        if u not in v_vertex_numbers:
            v_queue.append(reverse)

        self.adj[u] = u_queue
        self.adj[v] = v_queue
        self.num_edges += 1
        
        if cap > self.max_capacity:
            self.max_capacity = cap

    def get_capacities(self):
        capacities = np.zeros((self.num_vertices, self.num_vertices))
        for i in range(self.num_vertices):
            for j in range(self.num_vertices):
                if i != j:
                    try:
                        capacities[i, j] = self.edges[(i, j)].capacity
                    except KeyError:
                        pass
        return capacities

    def admissible(self, s, t):
        for i in range(len(self.level)):
            self.level[i] = -1
        self.level[s] = 0
        q = deque([])
        q.append(s)
        while len(q) > 0:
            u = q.popleft()
            for e in self.adj[u]:
                if self.level[e.v] < 0 and e.flow < e.capacity:
                    self.level[e.v] = self.level[u] + 1
                    q.append(e.v)
        return self.level[t] >= 0

    def push_flow(self, u, flow, t, start):
        if u == t:
            return flow
        initial = start[u]
        for i in range(initial, len(self.adj[u])):
            e = self.adj[u][start[u]]
            if self.level[e.v] == self.level[u] + 1 and e.flow < e.capacity:
                current_node = min(flow, e.capacity - e.flow)
                temp = self.push_flow(e.v, current_node, t, start)
                if temp > 0:
                    e.flow += temp
                    self.adj[e.v][e.reverse].flow -= temp
                    return temp
            start[u] += 1
        return 0

    def dinic(self, s, t):
        if s == t:
            return -1
        total = 0
        while self.admissible(s, t):
            start = [0] * (self.num_vertices + 1)
            flow = self.push_flow(s, sys.maxsize, t, start)
            while flow:
                total += flow
                flow = self.push_flow(s, sys.maxsize, t, start)
        return total

    # Goldberg-Rao binary length function
    def binary_length_function(self, i, j, Delta):
        if self.resid_cap[i, j] >= 3 * Delta:
            return 0
        else:
            return 1

    # Get distance labels for all vertices
    def get_distance_labels(self, s, t, Delta, distance_labels):
        visited = np.zeros((self.num_vertices, self.num_vertices))
        distance_labels = [float('inf')] * self.num_vertices

        # Use DFS to find distance from each vertex to sink with
        # respect to the binary length function
        def DFS_for_distances(v):
            if v == t:
                distance_labels[t] = 0
            elif distance_labels[v] == float('inf'):
                for u in self.adj[v]:
                    u = u.v
                    if self.resid_cap[v, u] > 0 and visited[u, v] == 0:
                        visited[u, v] = 1
                        bin_length = self.binary_length_function(v, u, Delta)
                        distance_labels[v] = min(DFS_for_distances(u) + bin_length, distance_labels[v])
            return distance_labels[v]

        DFS_for_distances(s)
        return distance_labels
    
    def blocking_flow(self, s, t, distance_labels, Delta, flow_matrix):
        visited = np.zeros((self.num_vertices, 1))
        source_vertex = [-420] * self.num_vertices
        mincap = [float('inf')] * self.num_vertices
        
        node_stack = [s]
        while node_stack:
            u = node_stack.pop()
            visited[u, 0] = 1
            if u == t:
                break
            for v in self.adj[u]:
                v = v.v
                bin_length = self.binary_length_function(u, v, Delta)
                admissible = distance_labels[u] == distance_labels[v] + bin_length
                # Construct admissible graph and contract SCCs induced by zero-length arcs
                if not visited[v, 0] and self.resid_cap[u, v] > 0 and admissible:
                    node_stack.append(v)
                    source_vertex[v] = u
                    mincap[v] = min(mincap[u], self.resid_cap[u, v])

        # Push flow
        if mincap[t] != float('inf'):
            current_node = t
            while current_node:
                self.resid_cap[source_vertex[current_node], current_node] -= mincap[t]
                flow_matrix[source_vertex[current_node], current_node] += mincap[t]
                self.resid_cap[current_node, source_vertex[current_node]] += mincap[t]
                flow_matrix[current_node, source_vertex[current_node]] -= mincap[t]
                current_node = source_vertex[current_node]
        return flow_matrix

    def rao(self, s, t):
        self.resid_cap = self.get_capacities()
        flow_matrix = np.zeros((self.num_vertices, self.num_vertices), dtype=np.int)
        # F = mU is initially an upper bound on total flow
        F = self.num_edges * self.max_capacity
        Lambda = min(self.num_edges ** .5, self.num_vertices ** 1.5)
        distance_labels = [float('inf')] * self.num_vertices
        while F >= 1:
            Delta = math.ceil(F / (2 * Lambda))
            for count in range(1, math.ceil(5 * Lambda)):
                # Get distance labelings
                distance_labels = self.get_distance_labels(s, t, Delta, distance_labels)
                # Construct admissible graph, contract SCCs induced by 
                # zero-length arcs, and run blocking flow.
                flow_matrix = self.blocking_flow(0, 5, distance_labels, Delta, flow_matrix)
            F = F/2
        print("Rao-Goldberg flow matrix:")
        print(flow_matrix)
        return sum(flow_matrix[:, t])

print("Graph 1")
general = Graph(6, [])
general.add_edge(0, 1, 16)
general.add_edge(0, 2, 13)
general.add_edge(1, 2, 10)
general.add_edge(1, 3, 12)
general.add_edge(2, 1, 4)
general.add_edge(2, 4, 14)
general.add_edge(3, 2, 9)
general.add_edge(3, 5, 20)
general.add_edge(4, 3, 7)
general.add_edge(4, 5, 4)

blocking_flow = general.dinic(0, 5)
print("Dinic blocking max flow =", blocking_flow)
print("Rao-Goldberg max flow =", general.rao(0, 5))

print("\nGraph 2")
general = Graph(6, [])
general.add_edge(0, 1, 10)
general.add_edge(0, 2, 10)
general.add_edge(1, 2, 2)
general.add_edge(1, 3, 4)
general.add_edge(1, 4, 8)
general.add_edge(2, 4, 9)
general.add_edge(3, 5, 10)
general.add_edge(4, 3, 6)
general.add_edge(4, 5, 10)

blocking_flow = general.dinic(0, 5)
print("Dinic blocking max flow =", blocking_flow)
print("Rao-Goldberg max flow =", general.rao(0, 5))
