# Adelaide CBD - Chinese Postman Route Generator
# This script calculates the shortest possible walking route that traverses every street
# in the Adelaide CBD (North, West, South, and East Terraces) at least once.
# The solution is visualized on an interactive HTML map.

import osmnx as ox
import networkx as nx
import folium
from itertools import combinations
import time

def get_adelaide_cbd_graph():
    """
    Downloads and prepares the street network graph for Adelaide's CBD.
    The bounding box is defined by the four main terraces.
    """
    print("Downloading street map for Adelaide CBD from OpenStreetMap...")
    # Bounding box for Adelaide CBD (North, South, East, West Terraces)
    # A slightly expanded bounding box can help prevent "no nodes found" errors
    # by ensuring all intersections on the terraces are included.
    # North, South, East, West
    north, south, east, west = -34.9170, -34.9360, 138.6115, 138.5870
    
    # Download the street network data
    # Using 'walk' network type is most appropriate for this problem.
    # The osmnx API now expects the bounding box coordinates as a single tuple.
    graph = ox.graph_from_bbox((north, south, east, west), network_type='walk')
    
    # Project the graph to a local UTM zone to get distances in meters
    graph_proj = ox.project_graph(graph)
    
    print("Street network downloaded and projected successfully.")
    return graph_proj

def solve_chinese_postman(graph):
    """
    Solves the Chinese Postman Problem for the given graph.
    This involves finding nodes with an odd number of connections,
    finding the shortest paths between them, and creating an optimal
    set of "duplicate" paths to make the graph Eulerian.
    """
    print("Solving the Chinese Postman Problem...")
    
    # 1. Identify nodes with an odd degree (odd number of streets connected)
    odd_degree_nodes = [node for node, degree in graph.degree() if degree % 2 != 0]
    print(f"Found {len(odd_degree_nodes)} odd-degree nodes. These intersections will need to be re-visited.")

    # 2. Compute all-pairs shortest paths between the odd-degree nodes
    print("Calculating shortest paths between all pairs of odd-degree nodes...")
    start_time = time.time()
    # This creates pairs of odd nodes, e.g., (A, B), (A, C), (B, C)
    odd_node_pairs = list(combinations(odd_degree_nodes, 2))
    
    # Use a dictionary to store the costs (path lengths)
    costs = {}
    for pair in odd_node_pairs:
        # We use Dijkstra's algorithm to find the shortest path length.
        # The 'weight' attribute is the length of the street segment in meters.
        costs[pair] = nx.dijkstra_path_length(graph, pair[0], pair[1], weight='length')
    
    print(f"Path calculations finished in {time.time() - start_time:.2f} seconds.")

    # 3. Find a 'perfect matching' with minimum weight
    # This is the most crucial step: we want to pair up the odd nodes in a way
    # that the total distance of the paths connecting them is minimized.
    # We create a new complete graph where nodes are the odd nodes and edges
    # are weighted by the shortest path distance.
    print("Finding the optimal pairing of odd-degree nodes...")
    g_odd_complete = nx.Graph()
    g_odd_complete.add_weighted_edges_from([(pair[0], pair[1], cost) for pair, cost in costs.items()])
    
    # We use a max_weight_matching algorithm. To find the minimum, we invert the weights.
    # A higher negative number is "smaller", so the algorithm will find the minimum sum.
    for u, v, data in g_odd_complete.edges(data=True):
        data['weight'] = -data['weight']

    # This algorithm finds the pairing that maximizes the (negative) weight,
    # which is equivalent to minimizing the original positive path lengths.
    min_weight_matching = nx.max_weight_matching(g_odd_complete, maxcardinality=True)
    print("Optimal pairing found.")

    # 4. Augment the original graph
    # We add the shortest paths from the matching to the original graph to create
    # a new graph that has an Eulerian circuit. We use a MultiGraph to allow
    # for parallel edges (representing streets walked more than once).
    print("Creating the final augmented graph for the route...")
    graph_augmented = nx.MultiGraph(graph.copy())
    
    total_added_distance = 0
    for pair in min_weight_matching:
        path = nx.dijkstra_path(graph, pair[0], pair[1], weight='length')
        total_added_distance += costs[(min(pair), max(pair))] # Add to running total
        
        # Add the edges of this shortest path to the augmented graph
        for i in range(len(path) - 1):
            u, v = path[i], path[i+1]
            # It's important to get the original edge data (like length)
            edge_data = graph_augmented.get_edge_data(u, v)
            if edge_data:
                 # If multiple edges exist, pick the first one's data
                first_key = list(edge_data.keys())[0]
                graph_augmented.add_edge(u, v, **edge_data[first_key])
            else:
                 # This case should ideally not happen in a connected graph
                 print(f"Warning: Could not find edge data for {u}-{v}")


    # 5. Find the Eulerian circuit in the new augmented graph
    print("Calculating the final Eulerian circuit...")
    # An Eulerian circuit visits every edge exactly once. Because we added the
    # duplicated paths, this circuit now covers all original streets.
    start_node = list(graph_augmented.nodes())[0]
    eulerian_circuit = list(nx.eulerian_circuit(graph_augmented, source=start_node, keys=True))
    
    # Calculate total distance
    total_distance_orig = sum(data['length'] for u, v, data in graph.edges(data=True))
    total_distance_final = sum(data['length'] for u, v, data in graph_augmented.edges(data=True))

    print("\n--- Calculation Summary ---")
    print(f"Total length of all streets: {total_distance_orig / 1000:.2f} km")
    print(f"Extra distance to be re-walked: {total_added_distance / 1000:.2f} km")
    print(f"Total walk distance for the complete route: {total_distance_final / 1000:.2f} km")
    
    return eulerian_circuit, graph_augmented

def plot_on_map(graph, circuit, filename="adelaide_cbd_walk_route.html"):
    """
    Plots the calculated circuit on an interactive Folium map.
    """
    print(f"\nGenerating interactive map: {filename}...")
    
    # Get the non-projected version of the graph for lat/lon coordinates
    graph_unproj = ox.project_graph(graph, to_crs='epsg:4326')
    
    # Create a list of lat/lon points for the route
    route_coords = []
    for u, v, key in circuit:
        # Get the geometry (line) of the street segment
        edge_data = graph_unproj.get_edge_data(u, v, key)
        if 'geometry' in edge_data:
            # The geometry is a LineString, get its coordinates
            route_coords.extend(list(edge_data['geometry'].coords))
        else:
            # Fallback for edges without a geometry attribute
            u_node = graph_unproj.nodes[u]
            v_node = graph_unproj.nodes[v]
            route_coords.append((u_node['y'], u_node['x']))
            route_coords.append((v_node['y'], v_node['x']))

    # Create the map, centered on the route
    map_center = [route_coords[0][0], route_coords[0][1]]
    m = folium.Map(location=map_center, zoom_start=15, tiles='cartodbpositron')

    # Add the route as a PolyLine
    folium.PolyLine(
        route_coords,
        color='#e63946',  # A nice red color
        weight=3,
        opacity=0.8
    ).add_to(m)
    
    # Add start and end markers
    folium.Marker(
        location=route_coords[0],
        popup='Start/End Point',
        icon=folium.Icon(color='green', icon='play')
    ).add_to(m)

    # Save the map to an HTML file
    m.save(filename)
    print("Map generated successfully!")


if __name__ == "__main__":
    # Main execution block
    cbd_graph = get_adelaide_cbd_graph()
    
    # Check if the graph is connected, which is a requirement
    if not nx.is_connected(cbd_graph):
        print("Error: The graph is not connected. Cannot solve.")
    else:
        circuit, augmented_graph = solve_chinese_postman(cbd_graph)
        plot_on_map(augmented_graph, circuit)
        print(f"\nProcess complete. Open 'adelaide_cbd_walk_route.html' in a web browser to see the route.")


