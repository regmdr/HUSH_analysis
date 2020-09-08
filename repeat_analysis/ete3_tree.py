import numpy as np  
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace

def layout(node):
    if node.is_leaf():
        # Add node name to laef nodes
        N = AttrFace("name", fsize=11, fgcolor="black")
        faces.add_face_to_node(N, node, 0)
    if "weight" in node.features:
        # Creates a sphere face whose size is proportional to node's
        # feature "weight"
        C = CircleFace(radius=node.weight, color="blue", style="circle")
        # Let's make the sphere transparent
        C.opacity = 0.3
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="branch-top") #position= "float-behind") #position="float")


def read_weights(file):
    weights ={}
    with open(file, "r") as f:
        for line in f:
            values = line.strip().split(",")
            weights[values[0]] = values[1]
    return weights



############################################################

tree_weights =read_weights("tree_weights.proportion.txt")



test = Tree(newick="curated_tree.txt")
for n in test.traverse():
    if n.name in tree_weights:
        n.add_features(weight=np.float(tree_weights[n.name]))
    else:
        n.add_features(weight=np.float(0))


# Checking order
for n in test.traverse():
	print(n.weight, n.name)# .get_leaf_names()


# Create an empty TreeStyle
ts = TreeStyle()

# Set our custom layout function
ts.layout_fn = layout

# Draw a tree
ts.mode = "r"

# We will add node names manually
ts.show_leaf_name = False
# Show branch data
ts.show_branch_length = False
ts.show_branch_support = False
ts.optimal_scale_level=True    
ts.aligned_table_style=True   
# ts.complete_branch_lines_when_necessary = True

      

# test.render("circle_map.genomic.pdf", w=10000, dpi=1000, tree_style=ts)
# test.render("circle_map.up_families.M-P.pdf", w=10000, dpi=1000, tree_style=ts)
test.render("circle_map.pdf", w=10000, dpi=1000, tree_style=ts)
