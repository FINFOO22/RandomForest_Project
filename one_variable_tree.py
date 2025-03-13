#!/usr/bin/python3

# Install Packages
from numpy.random import seed
from numpy.random import normal
from matplotlib.pyplot import hist
from matplotlib.pyplot import scatter
import matplotlib.pyplot as plt
from one_variable_tree_functions import *





dependent_variable, independent_variable = stable_normally_distributed_plot(100, 20)
#print(dependent_variable)
#print(independent_variable)

average, tree = first_ssr_calculation(dependent_variable=dependent_variable, independent_variable=independent_variable)
#print(average)
#print(tree)


# I need code that acknowledges where I am in the tree AND if it is a dead end
# for example 0 is left, 1 is right
# dead end = FALSE
# code will always attempt to go left until dead end = TRUE
# then will go right THEN see if it can go left again 
# if it is Right (1) and dead end = TRUE
# it will go back up and try right , then left directly after
# make a test tree
test_tree = [[[5, 4, 2, 5, 3, 4],[4, 5, 6, 7, 8, 9]],
             {5.3:
              [{2.1:[[1],[10]]},
               [[17,18,23,19,19],[11,12,13,14,15]]]}]
test_tree_actual = tree

pathway_list, dead_end, route_dictionary, what_is_it = navigate_tree(test_tree_actual, pathway_list=[1])
print(pathway_list)
print(dead_end)
print(route_dictionary)
print(what_is_it)




# 
# for example, it will now know that [1,0,1] is a dead end and proceed to the next theoretical option which would be [1,1] (as it is left to right, [1,0,0] wouldve already been done)
# make dictionaru   
dead_end_dict = {'dead': {}, 'open': {'top': average}} # key can be a tuple of the route and value is what it is
route = [0]
while True:
    pathway_list, dead_end, route_dictionary, what_is_it = navigate_tree(test_tree_actual, pathway_list=route)
    if what_is_it == 'branch':
        print(f'route {route} still requires completing')
        # run a function that looks at these data points and completes what it is meant to be
        # so if it is less than 20/7 data points, get average of the independent variable and thats the dead end leaf
        # IF not, then perform the ssr again
        print(route_dictionary[tuple(route)])
        if len(route_dictionary[tuple(route)][0]) < 7:
            leaf_value = sum(route_dictionary[tuple(route)][0])/len(route_dictionary[tuple(route)][0])
            leaf_dict = dead_end_dict['dead'] 
            leaf_dict[tuple(route)] = leaf_value
            dead_end_dict['dead'] = leaf_dict
    elif what_is_it == 'dict':
        print(f'route {route} is the start of another branch')
        # add this to the open key of dead_end_dict


    break
print(dead_end_dict)

# next - script looks at route, sees if it is in the dead, or open keys of the dictionary 
# if in the dead move on to the next suitable path , so if [0,0], change it to [0,1]
# then look to see if that is in open or closed
# if in open, use the navigate_tree function to obtain what it is then perform the next task whether it is a dictionary,
# or a branch etc 


#









quit()






[[[],[]],[[],[]]]