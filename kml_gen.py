# street_network_to_kml.py
#
# Description:
# This script fetches street network data (roads, paths, etc.) from the
# OpenStreetMap Overpass API for a specified latitude/longitude bounding box.
# It then processes this data and generates a KML file, where each street is
# represented as a LineString. This KML file can be opened in applications like
# Google Earth Pro or Google My Maps for visualization and route planning.
#
# Prerequisites:
# - Python 3
# - 'requests' library: pip install requests
# - 'simplekml' library: pip install simplekml

import requests
import simplekml
import sys

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
    # The Overpass QL query to retrieve all ways designated as 'highway'
    # and all the nodes that make up those ways.
    bbox_str = f"{bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}"
    overpass_query = f"""
    [out:json][timeout:30];
    (
      way({bbox_str})[highway][name];
    );
    (._;>;);
    out;
    """

    print("Sending request to Overpass API. This may take a moment for larger areas...")
    try:
        response = requests.post(OVERPASS_URL, data={'data': overpass_query})
        # Raise an HTTPError for bad responses (4xx or 5xx)
        response.raise_for_status()
        print("Data received successfully.")
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from Overpass API: {e}", file=sys.stderr)
        return None

def create_kml_file(data, output_filename="street_network.kml"):
    """
    Creates a KML file from the structured OpenStreetMap data.

    Args:
        data (dict): The OSM data retrieved from the Overpass API.
        output_filename (str): The name of the KML file to be saved.
    """
    if not data or 'elements' not in data:
        print("No data available to create KML file.", file=sys.stderr)
        return

    print(f"Processing data and generating {output_filename}...")

    # A dictionary to hold node IDs and their coordinates for quick access.
    nodes = {
        element['id']: (element['lon'], element['lat'])
        for element in data['elements'] if element['type'] == 'node'
    }

    kml = simplekml.Kml(name="Street Network", description="All streets within the specified bounding box.")

    # Iterate through all 'way' elements, which represent streets.
    for element in data['elements']:
        if element['type'] == 'way':
            tags = element.get('tags', {})
            street_name = tags.get('name', 'Unnamed Street')
            highway_type = tags.get('highway', 'unknown')

            # Collect the coordinates for each node in the way.
            coords = [nodes[node_id] for node_id in element.get('nodes', []) if node_id in nodes]

            # A line needs at least two points to be valid.
            if len(coords) >= 2:
                linestring = kml.newlinestring(name=street_name)
                linestring.coords = coords
                linestring.description = f"Street Type: {highway_type}"

                # Apply some basic styling.
                linestring.style.linestyle.color = simplekml.Color.blue
                linestring.style.linestyle.width = 3

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
    #
    # You can get coordinates from a site like https://www.openstreetmap.org/
    # by navigating to your area, clicking 'Export', and manually selecting a region.
    # The coordinates will be displayed on the left panel.
    #
    # NOTE: Keep the area reasonably small (e.g., a neighborhood or small town)
    # to avoid overloading the free Overpass API.
    # Example for a small area in Adelaide, Australia.
    bounding_box = (-34.93702364775846, 138.58745115505553, -34.91970776999262, 138.61716446913886)
    # -34.91970776999262, 138.61716446913886
    # -34.93702364775846, 138.58745115505553

    print("-" * 50)
    print("OpenStreetMap Street Network to KML Generator")
    print("-" * 50)
    print(f"Using Bounding Box (S, W, N, E): {bounding_box}")

    osm_data = fetch_street_data(bounding_box)

    if osm_data:
        create_kml_file(osm_data)

if __name__ == "__main__":
    main()
