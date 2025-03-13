#!/usr/bin/python3

# Install Packages
from numpy.random import seed
from numpy.random import normal
from matplotlib.pyplot import hist
from matplotlib.pyplot import scatter
import matplotlib.pyplot as plt
import random
from random import randint

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
        if i > 0.05 and i < 0.995:
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


def bootstrap(data, no_times=1):
    #
    # print(data)
    no_trees = 20
    # this would be done usually via obtaining the uploaded datas keys OR column headings if the data s in pandas format
    bootstrapped_data = {'bioavailability': None,
    'other_biomarker': [],
    'treated': [],
    'mgkg': [],
    'biomarker': []}
    indexes = [randint(0, (len(data[list(data.keys())[0]])-1)) for i in range(0, len(data[list(data.keys())[0]]))]
    unused_samples = [i if i not in indexes else None for i in range(0, len(data[list(data.keys())[0]]))]
    unused_samples = [x for x in unused_samples if x is not None]
    #print(unused_samples)
    #print(indexes)

    # here i have obtained a system that ranomdly selects, and NOTES which indexes werent used in the Bootstrap
    # Now, reconstruct the dictionary (bootstrap_data)
    #print(bootstrapped_data)
    for i in list(bootstrapped_data.keys()):
        value = data[i]
        list_ = []
        for x in indexes:
            real = value[x]
            list_ = list_ + [data[i][x]]
        bootstrapped_data[i] = list_

    #print(bootstrapped_data)
    return bootstrapped_data





def best_variable_ssr_calculation(dependent_variable, independent_variable, rearranged_indexes):
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
    ssr_ = ssr_list[index]
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
    split_rearranged_indexes = [[x for x in list(rearranged_indexes[:index+1])], [x for x in list(rearranged_indexes[index+1:])]]
    # Currently this code can produce the first split.
    return average, tree, ssr_, split_rearranged_indexes


## Now, have to order lists ascending to measure the ssrs (dependent variable in ascending order and the)
# compare the ssr for each variable.
# We will assume that dependent (biomarker) variable is the one that we will measure and look at sum squared residuals for each one and see which is the lowest
# practise with mgkg
def re_order(bootstrapped, ind_name, dep_name):
    # ordered_data = bootstrapped
    dep_var = bootstrapped[dep_name] # ACKNOWLEDGE this needs to change for when practise data no longer used
    ind_var = bootstrapped[ind_name] # need to order values based on the ind variable
    indexes = [x for x in range(0, len(ind_var))]
    rearranged_indexes = []
    for i in range(0,len(indexes)):
        try:
            if i == 0:
                rearranged_indexes = rearranged_indexes + [indexes[i]]
                continue
            elif ind_var[i] > ind_var[rearranged_indexes[-1]]:
                rearranged_indexes = rearranged_indexes + [i]
                continue
        except:
            pass
        for ind in range(0,len(rearranged_indexes)):
            if ind_var[i] <= ind_var[rearranged_indexes[ind]]:
                split_1 = rearranged_indexes[0:ind]
                split_2 = rearranged_indexes[ind:]
                rearranged_indexes = split_1 + [i] + split_2
                break
            elif ind_var[i] > ind_var[rearranged_indexes[ind]] and ind == (len(rearranged_indexes) - 1):
                rearranged_indexes = rearranged_indexes + [i]
    # print(rearranged_indexes)
    reorder_ind = [ind_var[x] for x in rearranged_indexes]
    #print(ind_name)
    #print(f'INDEPENDENT VAR {reorder_ind}')
    reorder_dep = [dep_var[x] for x in rearranged_indexes]
    #print(dep_name)
    #print(f'DEPENDENT VAR {reorder_dep}')
    return rearranged_indexes, reorder_dep, reorder_ind

def best_ind_var(bootstrapped, dep_name):
    names = list(bootstrapped.keys())
    root_choice = {}
    for n in names:
        if n == dep_name:
            pass
        else:
            rearranged_indexes, reorder_dep, reorder_ind = re_order(bootstrapped=bootstrapped, ind_name=n, dep_name=dep_name)
            average, tree, ssr_, split_rearranged_indexes = best_variable_ssr_calculation(dependent_variable=reorder_dep, 
                                                                independent_variable=reorder_ind,
                                                                rearranged_indexes=rearranged_indexes)
            root_choice[n] = [tree, split_rearranged_indexes, ssr_]
    return root_choice


def initial_multi_root_dict(root_choice, bootstrapped):
    ssr_ = []
    for key in list(root_choice.keys()):
        if len(ssr_) == 0:
            ssr_ = ssr_ + [key]
            ssr_ = ssr_ + [root_choice[key][2]]
        elif ssr_[1] > root_choice[key][2]:
            #print(key)#
            #print(root_choice[key][2])
            ssr_[0] = key
            ssr_[1] = root_choice[key][2]
    #print(ssr_)

    chosen_indexes = root_choice[ssr_[0]][1]
    rearranged_dict_of_vars = {}
    multi_tree = {'dead': {}, 'open': {'top': [ssr_[0], list(root_choice[ssr_[0]][0].keys())[0]]}}
    open_side = multi_tree['open']
    for key in list(bootstrapped.keys()):#
        var_ = bootstrapped[key]
        left = [var_[x] for x in chosen_indexes[0]]
        right = [var_[x] for x in chosen_indexes[1]]
        # make left and right branch
        try:
            branch_dict = open_side[tuple([0])]
            branch_dict[key] = left
        except:    
            open_side[tuple([0])] = {key: left}
        try:
            branch_dict = open_side[tuple([1])]
            branch_dict[key] = right
        except:    
            open_side[tuple([1])] = {key: right}
    multi_tree['open'] = open_side
    return multi_tree
    print(multi_tree)


## THIS IS NOT COMPLETE YET, MORE TO DO 
def complete_MULTI_tree(root_dict, route=[0]):
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
        complete_MULTI_tree(root_dict=root_dict, route=route)
    
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
        
        complete_MULTI_tree(root_dict=root_dict, route=left)
        #complete_tree(root_dict=root_dict, route=left_route)

    return root_dict

