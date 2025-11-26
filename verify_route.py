import xml.etree.ElementTree as ET
import networkx as nx
import math
from solve_kml_cpp import parse_kml, build_graph, haversine_distance

# Constants
ORIGINAL_KML = r'd:\Development\cbd-map\cbd-map.kml'
ROUTE_KML = r'd:\Development\cbd-map\full_route.kml'

def parse_route_kml(file_path):
    """Parses the generated route KML to get the sequence of coordinates."""
    tree = ET.parse(file_path)
    root = tree.getroot()
    ns = {'kml': 'http://www.opengis.net/kml/2.2'}
    
    coords_text = root.find('.//kml:LineString/kml:coordinates', ns).text
    coords = []
    for coord in coords_text.strip().split():
        parts = coord.split(',')
        lon = float(parts[0])
        lat = float(parts[1])
        coords.append((lon, lat))
    return coords

def verify_coverage():
    print("Building original graph...")
    lines = parse_kml(ORIGINAL_KML)
    G = build_graph(lines)
    
    # We know the solver only works on the largest connected component
    if not nx.is_connected(G):
        largest_cc = max(nx.connected_components(G), key=len)
        G_component = G.subgraph(largest_cc).copy()
        print(f"Original graph has {len(G)} nodes. Verifying against largest component with {len(G_component)} nodes.")
    else:
        G_component = G
        
    print("Parsing generated route...")
    route_coords = parse_route_kml(ROUTE_KML)
    print(f"Route has {len(route_coords)} points.")
    
    # Check coverage
    # We need to see if every edge in G_component is traversed by route_coords.
    # This is tricky because of floating point precision.
    # We can check if every node in G_component is visited.
    # And we can check total length.
    
    # 1. Node coverage
    # Snap route coords to graph nodes
    def snap(coord):
        return (round(coord[0], 5), round(coord[1], 5))
    
    route_nodes = set(snap(c) for c in route_coords)
    graph_nodes = set(G_component.nodes())
    
    missing_nodes = graph_nodes - route_nodes
    print(f"Nodes in graph: {len(graph_nodes)}")
    print(f"Nodes visited by route: {len(route_nodes.intersection(graph_nodes))}")
    print(f"Missing nodes: {len(missing_nodes)}")
    
    if len(missing_nodes) > 0:
        print("WARNING: Some nodes were not visited!")
    else:
        print("SUCCESS: All nodes in the connected component were visited.")

    # 2. Length comparison
    total_graph_length = sum(d['weight'] for u, v, d in G_component.edges(data=True))
    
    route_length = 0
    for i in range(len(route_coords) - 1):
        route_length += haversine_distance(route_coords[i], route_coords[i+1])
        
    print(f"Total street length (component): {total_graph_length/1000:.2f} km")
    print(f"Total route length: {route_length/1000:.2f} km")
    
    if route_length >= total_graph_length:
        print("SUCCESS: Route length is sufficient to cover the graph.")
    else:
        print("WARNING: Route length is shorter than graph length (impossible for full coverage).")

if __name__ == "__main__":
    verify_coverage()
