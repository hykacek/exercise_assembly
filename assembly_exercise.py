def build_adjacency_list(edges, directed=False):
    graph = {}
    for u, v in edges:
        if u not in graph:
            graph[u] = []
        if v not in graph:
            graph[v] = []
    for u, v in edges:
        graph[u].append(v)
        if not directed:
            graph[v].append(u)
    return graph

def degrees(graph):
    degree_dict = {}
    for v in graph:
        degree_dict[v] = len(graph[v])
    return degree_dict

def in_out_degree(graph):
    in_degree = {v: 0 for v in graph}
    out_degree = {v: 0 for v in graph}
    for v in graph:
        out_degree[v] = len(graph[v])
        for neighbor in graph[v]:
            if neighbor not in in_degree:
                in_degree[neighbor] = 0
            in_degree[neighbor] += 1
    return in_degree, out_degree

edges = [("P1", "P2"), ("P1", "P3"), ("P2", "P3"), ("P3", "P4"), ("P4", "P5")]

def dfs(graph, start):
    visited = []
    stack = [start]
    while stack:
        v = stack.pop()
        if v not in visited:
            visited.append(v)
            for neighbor in sorted(graph.get(v, []), reverse=True):
                if neighbor not in visited:
                    stack.append(neighbor)
    return visited

def connected_components(graph):
    components = []
    visited = set()
    for v in sorted(graph.keys()):
        if v not in visited:
            component = dfs(graph,v)
            components.append(set(component))
            visited.update(component)
    return components

def is_weakly_connected(graph):
    shadow_edges = []
    for u, neighbors in graph.items():
        for v in neighbors:
            shadow_edges.append((u, v))
    shadow = build_adjacency_list(shadow_edges, directed=False)
    for v in graph:
        if v not in shadow:
            shadow[v] = []
    components = len(connected_components(shadow))
    return components <= 1

def eulerian_type(graph, directed=False):
    if not directed:
        deg = degrees(graph)
        non_zero = {v for v, d in deg.items() if d > 0}
        if non_zero:
            components = connected_components(graph)
            active_components = [c for c in components if c & non_zero]
            if len(active_components) > 1 or not non_zero.issubset(active_components[0]):
                return "none"

        odd_count = 0
        for v in deg:
            if deg[v] % 2 != 0:
                odd_count += 1
        
        if odd_count == 0:
            return "circuit"
        elif odd_count == 2:
            return "path"
        else:
            return "none"
            
    else:
        if not is_weakly_connected(graph):
            return "none"
            
        in_deg, out_deg = in_out_degree(graph)
        excess_plus = 0   
        excess_minus = 0  
        
        for v in graph:
            diff = out_deg[v] - in_deg[v]
            if diff == 0:
                continue
            elif diff == 1:
                excess_plus += 1
            elif diff == -1:
                excess_minus += 1
            else:
                return "none"
                
        if excess_plus == 0 and excess_minus == 0:
            return "circuit"
        elif excess_plus == 1 and excess_minus == 1:
            return "path"
        else:
            return "none"

def kmer_composition(text, k):
    kmers = []
    for i in range(len(text) - k + 1):
        kmers.append(text[i:i+k])
    return kmers

def debruijn_graph(kmers):
    edges = []
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        edges.append((prefix, suffix))

    return build_adjacency_list(edges, directed=True)
    
def spell_path(nodes):
    path = nodes[0]
    for node in nodes[1:]:
        path += node[-1]
    return path

def hierholzer(graph, start, verbose=False):
    unused_edges = {u: sorted(neighbors) for u, neighbors in graph.items()}
    walk = [start]
    
    while sum(len(edges) for edges in unused_edges.values()) > 0:
        for i, node in enumerate(walk):
            if unused_edges.get(node, []): 
                start_node = node
                sub_walk = [start_node]
                
                while unused_edges.get(start_node, []):
                    next_node = unused_edges[start_node].pop(0)
                    sub_walk.append(next_node)
                    start_node = next_node
                    
                walk = walk[:i] + sub_walk + walk[i+1:]
                break
                
    return walk

def eulerian_path(graph, verbose=False):
    start_node = None
    in_deg, out_deg = in_out_degree(graph)
    for v in graph:
        if out_deg[v] - in_deg[v] == 1:
            start_node = v
            break
        if start_node is None:
            start_node = min(graph)

    return hierholzer(graph, start_node, verbose=verbose)

def contigs(graph):
    in_deg, out_deg = in_out_degree(graph)
    contig_list = []
    not_1in_1out = []
    for v in graph:
        if in_deg[v] != 1 or out_deg[v] != 1:
            for neighbor in graph[v]:
                not_1in_1out.append((v, neighbor))
    for v, neighbor in not_1in_1out:
        walk = [v, neighbor]
        current = neighbor
        while in_deg[current] == 1 and out_deg[current] == 1:
            next_node = graph[current][0]
            walk.append(next_node)
            current = next_node
        contig_list.append(spell_path(walk))
        