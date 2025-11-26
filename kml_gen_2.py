# street_network_to_kml.py
#
# Description:
# This script fetches street network data from the OpenStreetMap Overpass API
# for a specified latitude/longitude bounding box. It intelligently processes
# the data to identify all intersections and splits the streets into individual
# segments between those intersections. The resulting KML is organized into
# folders for each street name, providing a clean and topologically correct
# representation of the street network graph for route planning.
#
# Prerequisites:
# - Python 3
# - 'requests' library: pip install requests
# - 'simplekml' library: pip install simplekml

import requests
import simplekml
import sys
from collections import defaultdict

# The public Overpass API endpoint.
OVERPASS_URL = "https://overpass-api.de/api/interpreter"

def fetch_street_data(bbox):
    """
    Fetches street data from the Overpass API for a given bounding box.

    Args:
        bbox (tuple): A tuple containing the bounding box coordinates in the
                      format (south_lat, west_lon, north_lat, east_lon).

    Returns:
        dict: A dictionary containing the parsed JSON response from the API,
              or None if an error occurred.
    """
    bbox_str = f"{bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}"
    # This query finds ways with a 'highway' tag, excluding common non-drivable paths.
    # The '!~' operator means 'does not match regex'.
    overpass_query = f"""
    [out:json][timeout:60];
    (
      way["highway"!~"^(footway|path|cycleway|steps|pedestrian|track|service)$"]({bbox_str});
    );
    (._;>;);
    out;
    """

    print("Sending request to Overpass API. This may take a moment for larger areas...")
    try:
        response = requests.post(OVERPASS_URL, data={'data': overpass_query})
        response.raise_for_status()
        print("Data received successfully.")
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from Overpass API: {e}", file=sys.stderr)
        return None

def create_kml_file(data, output_filename="street_network.kml"):
    """
    Creates a KML file from the structured OpenStreetMap data, breaking ways
    into segments at each intersection.

    Args:
        data (dict): The OSM data retrieved from the Overpass API.
        output_filename (str): The name of the KML file to be saved.
    """
    if not data or 'elements' not in data:
        print("No data available to create KML file.", file=sys.stderr)
        return

    print(f"Processing data and generating {output_filename}...")

    # Step 1: Pre-process nodes and ways
    nodes = {el['id']: (el['lon'], el['lat']) for el in data['elements'] if el['type'] == 'node'}
    ways = [el for el in data['elements'] if el['type'] == 'way']

    # Step 2: Identify all intersection nodes
    node_usage_count = defaultdict(int)
    for way in ways:
        for node_id in way.get('nodes', []):
            node_usage_count[node_id] += 1

    intersection_node_ids = set()
    for node_id, count in node_usage_count.items():
        if count > 1:
            intersection_node_ids.add(node_id)

    # Also, treat the start and end points of any way as "intersections"
    # to ensure all segments are correctly terminated.
    for way in ways:
        way_nodes = way.get('nodes', [])
        if len(way_nodes) > 0:
            intersection_node_ids.add(way_nodes[0])
            intersection_node_ids.add(way_nodes[-1])

    # Step 3: Build KML, segmenting ways at each intersection
    kml = simplekml.Kml(name="Street Network", description="All streets, segmented by intersection.")
    street_folders = {}

    for way in ways:
        tags = way.get('tags', {})
        street_name = tags.get('name', 'Unnamed Street')
        highway_type = tags.get('highway', 'unknown')
        
        # Create a new folder for the street if it doesn't exist
        if street_name not in street_folders:
            street_folders[street_name] = kml.newfolder(name=street_name)
        folder = street_folders[street_name]

        way_node_ids = way.get('nodes', [])
        if len(way_node_ids) < 2:
            continue

        segment_start_index = 0
        for i in range(1, len(way_node_ids)):
            node_id = way_node_ids[i]
            
            # A segment ends if the current node is an intersection,
            # or it's the very last node of the way.
            is_intersection = node_id in intersection_node_ids
            is_end_of_way = (i == len(way_node_ids) - 1)

            if is_intersection or is_end_of_way:
                # Extract the node IDs for this specific segment
                segment_node_ids = way_node_ids[segment_start_index : i + 1]
                
                if len(segment_node_ids) > 1:
                    coords = [nodes[nid] for nid in segment_node_ids if nid in nodes]
                    
                    if len(coords) > 1:
                        # Create a new linestring for this segment
                        linestring = folder.newlinestring(name=f"{street_name} segment")
                        linestring.coords = coords
                        linestring.description = f"Street Type: {highway_type}"
                        linestring.style.linestyle.color = simplekml.Color.blue
                        linestring.style.linestyle.width = 3
                
                # The next segment will start from the current node's index
                segment_start_index = i

    try:
        kml.save(output_filename)
        print(f"\nSuccessfully saved KML file to: {output_filename}")
        print("You can now open this file in Google Earth or import it into Google My Maps.")
    except Exception as e:
        print(f"Error saving KML file: {e}", file=sys.stderr)

def main():
    """
    Main execution function.
    """
    # --- USER ACTION REQUIRED ---
    # Define your bounding box here.
    # The format is: (South Latitude, West Longitude, North Latitude, East Longitude)
    bounding_box = (-34.93, 138.59, -34.92, 138.60)

    print("-" * 50)
    print("OpenStreetMap Street Network to KML Generator")
    print("-" * 50)
    print(f"Using Bounding Box (S, W, N, E): {bounding_box}")

    osm_data = fetch_street_data(bounding_box)

    if osm_data:
        create_kml_file(osm_data)

if __name__ == "__main__":
    main()

