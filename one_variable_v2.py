#!/usr/bin/python3

# Install Packages
from numpy.random import seed
from numpy.random import normal
from matplotlib.pyplot import hist
from matplotlib.pyplot import scatter
import matplotlib.pyplot as plt
from one_variable_tree_functions import *


#dependent_variable, independent_variable = stable_normally_distributed_plot(100, 20)

#average, tree = first_ssr_calculation(dependent_variable=dependent_variable, independent_variable=independent_variable)

#pathway_list, dead_end, route_dictionary, what_is_it = navigate_tree(tree)


root_dict = initial_root_dict()

root_dict = complete_tree(root_dict)
print(root_dict)




## Next steps
#  introduce option to change the Min No. data points for another branch
#  Potentially make a drawing function to draw out the tree.
#  Start the multi independent variables tree.
