#!/usr/bin/python3

# Install Packages
from numpy.random import seed
from numpy.random import normal
from matplotlib.pyplot import hist
from matplotlib.pyplot import scatter
import matplotlib.pyplot as plt
import random

def stable_normally_distributed_plot(size_, no_bins):
    # Reproducible ND
    seed(1)
    # Here is a practise at constructing a one variable tree
    # I will use a simple, normally distributed set of data points
    data = normal(loc=0, scale=1, size=size_)
    count, bins, ignored = hist(data, no_bins)
    count = [float(x) for x in count]

    dep_vars = []
    dep_var = 1
    for i in count:
        dep_vars.append(dep_var)
        dep_var += 1
    # Here is the 'data' we have made
    plt.clf()
    #plot = scatter(dep_vars, count)
    #plt.show()
    #plt.clf()
    #print(dep_vars)
    #print(count)
    return dep_vars, count




def first_ssr_calculation(dependent_variable, independent_variable):
    tree = {}
    # Start of the SSR 
    # Find the smallest SUM of SQUARED RESIDUALS
    ssr_list = []
    for i in range(0, len(dependent_variable), 1):
        # split the root, then calculate the SSR for all the dat points
        split_1 = dependent_variable[:i+1]
        split_2 = dependent_variable[i+1:]
        # calculate ssr based on the average of either split
        squared_resids_1 = [(x - (sum(split_1)/len(split_1)))**2 for x in split_1]
        squared_resids_2 = [(x - (sum(split_2)/len(split_2)))**2 for x in split_2]

        # add ssr together
        ssr = sum(squared_resids_1) + sum(squared_resids_2)
        ssr_list.append(ssr)
    #print(f'here is the ssr list {ssr_list}')
    # obtain index for largest value in ssr_list
    for i in range(0,len(ssr_list),1):
        if ssr_list[i] == min(ssr_list):
            index = i
    #print(f'here is the index: {index}')
    #print(independent_variable)
    #print(independent_variable[index])
    #print(independent_variable[index+1])
    # obtain the independent variable value that the sum of squared residuals is lowest
    average = (independent_variable[index] + independent_variable[index+1])/2
    # going to make the dictionary (tree)
    #print(f'here is the average{average}')
    tree[average] = [[[float(x) for x in list(independent_variable[:index+1])], 
                      dependent_variable[:index+1]], 
                      [independent_variable[index+1:], dependent_variable[index+1:]]]
    # Currently this code can produce the first split.
    return average, tree



def navigate_tree(tree_start, pathway_list=[], dead_end=False, route_dictionary={}):
    if len(pathway_list) == 0:
        pathway_list.append(0)
    path = tree_start # this is the first two splits at the top of the tree
    pathway_list_progression = []

    for i in range(0,len(pathway_list),1):
        if isinstance(path, dict):
            path = path[list(path.keys())[0]]
        pathway_list_progression.append(pathway_list[i])

        path = path[pathway_list[i]] # left or right
        # based on the structure of the unfinished tree
        # this leaves a list of lists (the list of dependent values and a list of the independent values)
        # OR 
        # a dictionary
        # OR 
        # A list with a single value inside (this would be the FINAL leaf)
        if isinstance(path, dict):
            pass
        elif isinstance(path[0], float) or isinstance(path[0], int):
            dead_end = True
            # path = path[0]
            break
        elif isinstance(path[0], list):
            dead_end = True
            break

    # you have finished the pathway route
    # print(path)
    route_dictionary[tuple(pathway_list_progression)] = path
    if isinstance(path, dict):
        what_is_it = 'dict'
    elif isinstance(path[0], list):
        what_is_it = 'branch'
        dead_end = True
    elif isinstance(path[0], float) or isinstance(path[0], int):
        what_is_it = 'final_leaf'
        dead_end = True
        
    
    return pathway_list, dead_end, route_dictionary, what_is_it


def initial_root_dict(dependent_variable=None, independent_variable=None): #, initial_root_dict={'dead': {}, 'open': {}}):
    if dependent_variable is None and independent_variable is None:
        independent_variable, dependent_variable = stable_normally_distributed_plot(1000, 30)
    average, tree = first_ssr_calculation(dependent_variable=dependent_variable, independent_variable=independent_variable)
    root_dict = {'dead': {}, 'open': {'top': average}}
    #print(tree)
    routes = tree[average]
    #print(routes[0])
    #print(routes[1])
    root_dict_open = root_dict['open']
    root_dict_open[tuple([0])] = routes[0]
    root_dict_open[tuple([1])] = routes[1]
    root_dict['open'] = root_dict_open
    return root_dict



def complete_tree(root_dict, route=[0]):
    ## start at [0]
    root_dict_open = root_dict['open']
    #print(root_dict)
    print(f'HERE IS THE ROOT DICT SO FAR {root_dict}')
    route_value = root_dict_open[tuple(route)]
    if len(route_value[0]) < 7: # just using 7 currently as data size is small
        pathway_value = (sum(route_value[1])/len(route_value[1])) # average of the route
        dead = pathway_value
        root_dict['dead'][tuple(route)] = dead
        del root_dict['open'][tuple(route)]
        # now have to write code that switches to the next 'chronological' route path  
        # this is called recursion
        i = 1
        
        while True:
            #print(route[0:-i + 1])
            try:
                if route[-i] == 0: # if it is 0, the next option is to change the the same indexed item in the route list to 1
                    if i == 1:
                        route[-i] = 1
                        break
                    else:
                        #print(route)
                        route = route[0:-i+1]
                        #print(route)
                        route[-1] = 1
                        break
                else:
                    i += 1
            except IndexError:
                return root_dict
        complete_tree(root_dict=root_dict, route=route)
    
    elif len(route_value[0]) >= 7:
        # print(f'this {route_value} needs doing')
        # edit the ssr calculation function
        #print(f'suppsed dependent {route_value[0]}')
        #print(f'supposed independent {route_value[1]}')
        average_, tree_ = first_ssr_calculation(dependent_variable=route_value[1], independent_variable=route_value[0])
        #print(average_)
        left = route + [0]
        right = route + [1]

        #print(tree_)
        root_dict_open = root_dict['open']
        root_dict_open[tuple(route)] = average_
        root_dict_open[tuple(left)] = tree_[average_][0]
        root_dict_open[tuple(right)] = tree_[average_][1]
        root_dict['open'] = root_dict_open
        
        complete_tree(root_dict=root_dict, route=left)
        #complete_tree(root_dict=root_dict, route=left_route)

    return root_dict

def more_variables(no_samples):
    seed(1)
    random_ = [(random.random()/5) for i in range(0, no_samples)]
    values = [x for x in range(-10, (-10+no_samples))]
    values_adjusted = [(values[x]+random_[x]) for x in range(0, no_samples)]
    # print(values_adjusted)
    sigmoid_vals = sigmoid_funct(values_adjusted)
    # now make a classification
    effective = []

    for i in sigmoid_vals:
        #print(i)
        if i > 0.05 and i < 0.95:
            #print(i)
            if random.random() > 0.1:
                effective = effective + [1]
            else:
                effective = effective + [0]
        else:
            if random.random() < 0.1:
                effective = effective + [1]
            else:
                effective = effective + [0]
    random__ = [(random.random()) for i in range(0, no_samples)]

    return sigmoid_vals, random__, effective

                       
def sigmoid_funct(x_vals):
    sigmoid_vals = [(1/(1+(2**-x))) for x in x_vals]
    return sigmoid_vals