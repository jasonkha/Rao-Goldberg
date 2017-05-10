from collections import deque
import sys
import math
import numpy as np
import copy

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
        self.edges[(v, u)] = forward
        u_queue.append(forward)
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
                current = min(flow, e.capacity - e.flow)
                temp = self.push_flow(e.v, current, t, start)
                if temp > 0:
                    e.flow += temp
                    self.adj[e.v][e.reverse].flow -= temp
                    return temp
            start[u] += 1
        return 0

    def blocking_max_flow(self, s, t, old_graph):
        old = copy.deepcopy(old_graph)
        if s == t:
            return -1
        total = 0
        while self.admissible(s, t):
            start = [0] * (self.num_vertices + 1)
            flow = self.push_flow(s, sys.maxsize, t, start)
            while flow:
                total += flow
                flow = self.push_flow(s, sys.maxsize, t, start)
        return total, old

    def binary_length_function(self, i, j, residual_capacity, Delta):
        if residual_capacity[i, j] >= 3 * Delta:
            return 0
        else:
            return 1

    def distance(self, s, t, residual_capacity, Delta, distances, visited):
        if s == t:
            distances[t] = 0
        elif distances[s] == float('inf'):
            for u in self.adj[s]:
                a = u.v
                if residual_capacity[s, a] > 0 and visited[s, a] == 0:
                    visited[s, a] = 1
                    new_distance = (self.distance(a, t, residual_capacity, Delta, distances, visited)[a] 
                                   + self.binary_length_function(s, a, residual_capacity, Delta))
                    distances[s] = min(distances[s], new_distance)
        return distances

    def get_admissible_distances(self, s, t, residual_capacity, Delta):
        distances = [float('inf')] * self.num_vertices
        visited = np.zeros((self.num_vertices, self.num_vertices))
        distances = self.distance(s, t, residual_capacity, Delta, distances, visited)
        return distances

    def augment_flow_scc(self, s, t, residual_capacity, flow_matrix, arc, Delta):
        visited = [False for i in range(self.num_vertices)]
        source = [0 for i in range(self.num_vertices)]
        min_capacity = [float('inf') for _ in range(self.num_vertices)]
        dfs_stack = [s]
        while dfs_stack:
            v = dfs_stack.pop()
            visited[v] = True
            if v == t:
                break
            for u in self.adj[v]:
                a = u.v
                if (not visited[a] and residual_capacity[v, a] > 0 and arc[v] == arc[a] 
                    + self.binary_length_function(v, a, residual_capacity, Delta)):
                    dfs_stack.append(a)
                    source[a] = v
                    min_capacity[a] = min(min_capacity[v], residual_capacity[v, a])

        if min_capacity[t] != float('inf'):
            dest = t
            G_a = copy.deepcopy(self)
            while dest and G_a.admissible(s, t):
                flow, G_a = G_a.blocking_max_flow(s, t, G_a)
                dest = source[dest]
        return flow

    def rao(self, s, t, capacities):
        flow_matrix = np.zeros((self.num_vertices, self.num_vertices), dtype=np.int)
        F = self.num_edges * self.max_capacity
        flow = 0
        Lambda = min(self.num_edges ** .5, self.num_vertices ** 1.5)
        while F >= 1:
            Delta = math.ceil(F / 2 * Lambda)
            for count in range(1, math.ceil(5 * Lambda)):
                d_l = self.get_admissible_distances(s, t, capacities, Delta)
                flow = self.augment_flow_scc(s, t, capacities, flow_matrix, d_l, Delta)
            F = F/2
        return flow

general = Graph(6, []);
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
capacities = general.get_capacities()
blocking_flow, general = general.blocking_max_flow(0, 5, general)
print("Dinic blocking max flow =", blocking_flow)
print("Rao-Goldberg max flow =", general.rao(0, 5, capacities))
