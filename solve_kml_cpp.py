import xml.etree.ElementTree as ET
import networkx as nx
import math
from itertools import combinations
import shapely.geometry
import shapely.ops

def parse_kml(file_path):
    """Parses KML and returns a list of coordinate lists (lines)."""
    tree = ET.parse(file_path)
    root = tree.getroot()
    
    # Namespace handling
    ns = {'kml': 'http://www.opengis.net/kml/2.2'}
    
    lines = []
    
    for placemark in root.findall('.//kml:Placemark', ns):
        line_string = placemark.find('.//kml:LineString', ns)
        if line_string is not None:
            coordinates_text = line_string.find('kml:coordinates', ns).text
            if coordinates_text:
                coords = []
                for coord in coordinates_text.strip().split():
                    parts = coord.split(',')
                    lon = float(parts[0])
                    lat = float(parts[1])
                    coords.append((lon, lat))
                if len(coords) > 1:
                    lines.append(coords)
    return lines

def planarize_lines(lines):
    """
    Uses shapely to planarize the lines (split at intersections).
    Returns a list of coordinate lists.
    """
    print("Planarizing graph...")
    shapely_lines = [shapely.geometry.LineString(line) for line in lines]
    
    # unary_union merges overlapping lines and splits them at intersections
    merged = shapely.ops.unary_union(shapely_lines)
    
    planarized_lines = []
    if isinstance(merged, shapely.geometry.MultiLineString):
        for line in merged.geoms:
            planarized_lines.append(list(line.coords))
    elif isinstance(merged, shapely.geometry.LineString):
        planarized_lines.append(list(merged.coords))
    elif isinstance(merged, shapely.geometry.GeometryCollection):
         for geom in merged.geoms:
            if isinstance(geom, shapely.geometry.LineString):
                planarized_lines.append(list(geom.coords))
            elif isinstance(geom, shapely.geometry.MultiLineString):
                for line in geom.geoms:
                    planarized_lines.append(list(line.coords))
    
    print(f"Planarization complete. Converted {len(lines)} segments into {len(planarized_lines)} segments.")
    return planarized_lines

def haversine_distance(coord1, coord2):
    """Calculates Haversine distance between two (lon, lat) points in meters."""
    R = 6371000  # Earth radius in meters
    lon1, lat1 = coord1
    lon2, lat2 = coord2
    
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)
    
    a = math.sin(delta_phi / 2)**2 + \
        math.cos(phi1) * math.cos(phi2) * \
        math.sin(delta_lambda / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    return R * c

def build_graph(lines):
    """Builds a NetworkX MultiGraph from lines."""
    G = nx.MultiGraph()
    
    # Snap coordinates to a grid to handle precision issues
    # 5 decimal places is approx 1.1 meters
    def snap(coord):
        return (round(coord[0], 5), round(coord[1], 5))
    
    for line in lines:
        for i in range(len(line) - 1):
            u = snap(line[i])
            v = snap(line[i+1])
            
            if u != v:
                dist = haversine_distance(u, v)
                # Store original coords for reconstruction
                G.add_edge(u, v, weight=dist, coords=[line[i], line[i+1]])
                
    return G

def solve_cpp(G):
    """Solves the Chinese Postman Problem."""
    
    # 1. Check Connectivity
    if not nx.is_connected(G):
        print("Graph not connected. Using largest component.")
        largest_cc = max(nx.connected_components(G), key=len)
        G = G.subgraph(largest_cc).copy()
        print(f"Using component with {len(G)} nodes out of {len(G.nodes())}.")
    
    # 2. Identify Odd Degree Nodes
    odd_degree_nodes = [v for v, d in G.degree() if d % 2 == 1]
    print(f"Number of odd degree nodes: {len(odd_degree_nodes)}")
    
    # 3. Find Minimum Weight Matching
    # We need to add edges between odd nodes to make them even.
    # The weight of these edges is the shortest path distance.
    
    odd_node_pairs = list(combinations(odd_degree_nodes, 2))
    print(f"Calculating shortest paths for {len(odd_node_pairs)} pairs...")
    
    # Optimization: Use all_pairs_dijkstra only for odd nodes if possible, 
    # or just run it once if graph is small enough. 
    # For 400 nodes, all_pairs is fine.
    
    path_lengths = dict(nx.all_pairs_dijkstra_path_length(G, weight='weight'))
    
    G_odd = nx.Graph()
    for u, v in odd_node_pairs:
        if v in path_lengths[u]:
            G_odd.add_edge(u, v, weight=path_lengths[u][v])
            
    print("Calculating minimum weight matching...")
    matching = nx.min_weight_matching(G_odd, weight='weight')
    
    # 4. Augment the Graph
    print("Augmenting graph...")
    G_aug = G.copy()
    
    for u, v in matching:
        path = nx.dijkstra_path(G, u, v, weight='weight')
        # Add edges for the path
        for i in range(len(path) - 1):
            u_path = path[i]
            v_path = path[i+1]
            # Find the edge with the minimum weight between u_path and v_path
            # We need to copy the attributes
            edge_data = min(G[u_path][v_path].values(), key=lambda x: x['weight'])
            G_aug.add_edge(u_path, v_path, **edge_data)
            
    print(f"Graph augmented. Eulerian: {nx.is_eulerian(G_aug)}")
    
    # 5. Eulerian Circuit
    print("Calculating Eulerian circuit...")
    # keys=True is important to identify which parallel edge is used
    circuit = list(nx.eulerian_circuit(G_aug, source=list(G_aug.nodes())[0], keys=True))
    
    return circuit, G_aug

def create_kml_route(circuit, G, output_file):
    """Generates KML from the circuit."""
    kml = ET.Element('kml', xmlns="http://www.opengis.net/kml/2.2")
    doc = ET.SubElement(kml, 'Document')
    name = ET.SubElement(doc, 'name')
    name.text = "Adelaide CBD Walking Route"
    
    style = ET.SubElement(doc, 'Style', id="routeStyle")
    line_style = ET.SubElement(style, 'LineStyle')
    color = ET.SubElement(line_style, 'color')
    color.text = "ff0000ff" # Red
    width = ET.SubElement(line_style, 'width')
    width.text = "4"
    
    placemark = ET.SubElement(doc, 'Placemark')
    p_name = ET.SubElement(placemark, 'name')
    p_name.text = "Full Walking Route"
    p_style = ET.SubElement(placemark, 'styleUrl')
    p_style.text = "#routeStyle"
    
    line_string = ET.SubElement(placemark, 'LineString')
    tessellate = ET.SubElement(line_string, 'tessellate')
    tessellate.text = "1"
    coordinates = ET.SubElement(line_string, 'coordinates')
    
    # Reconstruct the full path coordinates
    full_coords_str = []
    
    for u, v, key in circuit:
        data = G[u][v][key]
        coords = data['coords']
        
        start_coord = coords[0]
        end_coord = coords[-1]
        
        dist_u_start = math.hypot(u[0]-start_coord[0], u[1]-start_coord[1])
        dist_u_end = math.hypot(u[0]-end_coord[0], u[1]-end_coord[1])
        
        if dist_u_end < dist_u_start:
            coords = coords[::-1]
            
        for lon, lat in coords:
            full_coords_str.append(f"{lon},{lat},0")
            
    coordinates.text = " ".join(full_coords_str)
    
    tree = ET.ElementTree(kml)
    ET.indent(tree, space="\t", level=0)
    tree.write(output_file, encoding='UTF-8', xml_declaration=True)
    print(f"KML route saved to {output_file}")

if __name__ == "__main__":
    input_kml = r'd:\Development\cbd-map\cbd-map.kml'
    output_kml = r'd:\Development\cbd-map\full_route.kml'
    
    print("Parsing KML...")
    lines = parse_kml(input_kml)
    print(f"Found {len(lines)} street segments.")
    
    # Planarize lines to fix connectivity
    lines = planarize_lines(lines)
    
    print("Building Graph...")
    G = build_graph(lines)
    print(f"Graph has {len(G)} nodes and {len(G.edges())} edges.")
    
    circuit, G_aug = solve_cpp(G)
    
    total_length = sum(G_aug[u][v][k]['weight'] for u, v, k in circuit)
    print(f"Route calculated with {len(circuit)} steps.")
    print(f"Total route distance: {total_length/1000:.2f} km")
    
    create_kml_route(circuit, G_aug, output_kml)
