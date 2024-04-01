import subprocess

def install_package(package_name):
    subprocess.check_call(['pip3', 'install', package_name])
packages=['ete3','pandas']
for package in packages:
    install_package(package)

from ete3 import Tree
import pandas as pd
def traverse_tree(tree, list_of_species):
    if len(tree.get_leaf_names()) < 2:
        return
    list_of_species.append(tree.get_leaf_names())
    for tr in tree.get_children():
        traverse_tree(tr, list_of_species)

def create_species_dataframe(lists_of_species):
    data = {}
    for i, species_list in enumerate(lists_of_species):
        for species in lists_of_species[0]:
            if species in species_list:
                if species not in data:
                    data[species] = [0] * len(lists_of_species)
                data[species][i] = 1

    df = pd.DataFrame(data, index=range(1, len(lists_of_species) + 1))
    df = df.T
    df = df.reset_index()
    df = df.rename(columns={'index': 'species_name'})
    return df
# Example Newick tree string
newick_str = input("Insert the species tree in newick format with an ending semi-colon: ")
tree = Tree(newick_str, format=1)
listf = []
traverse_tree(tree, listf)
df = create_species_dataframe(listf)
df.to_csv('results/tf-c-'+str(len(listf[0]))+'spcs.csv')

print("Saved tree_features in the parent directory")
