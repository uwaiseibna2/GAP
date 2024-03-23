import subprocess

def install_package(package_name):
    subprocess.check_call(['pip3', 'install', package_name])
packages=['ete3','pandas']
for package in packages:
    install_package(package)

from ete3 import Tree
import re
import pandas as pd
def generate_features(tree_str):
    tree = Tree(tree_str, format=1)

    # Assign names to nodes without names or with empty names
    for node in tree.traverse():
        if not node.name or node.name.strip() == "":
            node.name = f"Node_{id(node)}"

    # Initialize a dictionary to store features for each branch
    features = {}

    # Traverse the tree in postorder to ensure child branches are processed before their parent
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            features[node.name] = set()
            for descendant in node.iter_leaves():
                features[node.name].add(descendant.name)

    # Generate feature vectors
    species_set = {species for branch_features in features.values() for species in branch_features}
    feature_vectors = {}
    for branch, branch_features in features.items():
        feature_vector = [1 if species in branch_features else 0 for species in sorted(species_set)]
        feature_vectors[branch] = feature_vector

    return feature_vectors

# Example Newick tree string
newick_tree = input("Insert the species tree in newick format with an ending semi-colon: ")

features = generate_features(newick_tree)

species = re.findall(r'[a-zA-Z_ ]+', newick_tree)
tree_feats=pd.DataFrame(features)
tree_feats.index=species
tree_feats.columns=range(1, len(tree_feats.columns) + 1)
tree_feats.to_csv("tree_feats.csv")
print("Saved tree_features.csv in the parent directory")
