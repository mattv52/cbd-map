from fastkml import  kml
import sys
import networkx as nx

try:
    with open("test.kml", 'rt', encoding='utf-8') as f:
        doc = f.read().encode("utf-8")
    k = kml.KML()
    k = k.from_string(doc)
except Exception as e:
    print(f"Error reading or parsing KML file: {e}", file=sys.stderr)


# Create the KML object to store the parsed result
# Read in the KML string

# Next we perform some simple sanity checks

# Check that the number of features is correct
# This corresponds to the single ``Document``
graph = nx.MultiGraph()
# A dictionary to get all features from the KML, which are the LineStrings.
features = list(list(k.features())[0].features())
print(len(features))