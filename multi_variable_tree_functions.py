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
    try:
        average = (independent_variable[index] + independent_variable[index+1])/2
    except:
        average = independent_variable[index]
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
    sqrt = int((len(list(root_choice.keys()))**(1/2) + 1)) # THIS ROUNDS UP
    random_features = []
    i = 0 
    while i < sqrt:
        no = randint(0,sqrt)
        if no not in random_features:
            random_features = random_features + [no]
            i += 1
    random_features = [list(root_choice.keys())[x] for x in random_features]
    #print(random_features)
    #print(list(root_choice.keys()))
    #quit()
    for key in list(root_choice.keys()):
        if key in random_features:
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
def complete_MULTI_tree(root_dict, route=[0], dep_name='biomarker'):
    ## start at [0]
    root_dict_open = root_dict['open']
    #print(root_dict)
    print(f'HERE IS THE ROOT DICT SO FAR {root_dict}')
    route_value = root_dict_open[tuple(route)]
    # ABOVE give the dictionary of the variables
    # route_value = route_value[dep_name]

    if len(route_value[dep_name]) <= 7: # just using 7 currently as data size is small
        pathway_value = (sum(route_value[dep_name])/len(route_value[dep_name])) # average of the route
        dead = pathway_value
        root_dict['dead'][tuple(route)] = dead
        del root_dict['open'][tuple(route)]
        # now have to write code that switches to the next 'chronological' route path  
        print(f'Here we have UPDATED AND CLOSED A LEAF:')
        print(root_dict)
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
    
    elif len(route_value[dep_name]) > 7:
        # print(f'this {route_value} needs doing')
        # edit the ssr calculation function
        #print(f'suppsed dependent {route_value[0]}')
        #print(f'supposed independent {route_value[1]}')
        root_choice = best_ind_var(bootstrapped=route_value, dep_name=dep_name)
        #print(average_)


        print(f'HERE we have made it to splitting {route}, which leaves a root_choice of:')
        print(root_choice)
        #quit()
        continued_multi_root_dict(root_choice=root_choice, bootstrapped=route_value, root_dict=root_dict, route=route, dep_name='biomarker')

    return root_dict

def continued_multi_root_dict(root_choice, bootstrapped, root_dict, route, dep_name='biomarker'):
    ssr_ = []
    alternate_use = root_dict['open'][tuple(route)]
    route_ = route
    #print(root_choice['bioavailability'][2])
    sqrt = int((len(list(root_choice.keys()))**(1/2) + 1)) # THIS ROUNDS UP
    random_features = []
    i = 0 
    while i < sqrt:
        no = randint(0,sqrt)
        if no not in random_features:
            random_features = random_features + [no]
            i += 1
    random_features = [list(root_choice.keys())[x] for x in random_features]
    #quit()
    for key in list(root_choice.keys()):
        # select random features (sqrt of the total features)
        if key in random_features:
            try:
                if len(ssr_) == 0 and key in random_features:
                    ssr_ = ssr_ + [key]
                    ssr_ = ssr_ + [list(root_choice[key][0].keys())[0]]
                elif root_choice[ssr_[0]][2] >= root_choice[key][2] and key in random_features:
                    #print(key)#
                    #print(root_choice[key][2])
                    ssr_[0] = key
                    ssr_[1] = list(root_choice[key][0].keys())[0]
            except:
                print(root_choice[ssr_[0]][2])
                print(list(root_choice[key][0].keys()))
                print('hello') ## there is something wrong here

                quit()
        #print(ssr_[1])
        #quit()
    
    #print(ssr_)
    ssr_ = ssr_ + [root_choice[ssr_[0]][0][ssr_[1]][0][0], root_choice[ssr_[0]][0][ssr_[1]][1][0]] + root_choice[ssr_[0]][1]
    print(root_choice)
    print(ssr_)
    print(list(root_choice[ssr_[0]][0].values())[0][0][1])
    print(list(root_choice[ssr_[0]][0].values())[0][1][1])
    #quit()
    #quit()
    #print(ssr_)
    #quit()
    # here , calculate variance reduction FOR this split
    # variance reduction score = (variance_combined - (variance_left.proportion_left + variance_right.proportion_right)/variance_combined)
    #CALCULATE VARIANCE
    
    sample_left = list(root_choice[ssr_[0]][0].values())[0][0][1]
    sample_right = list(root_choice[ssr_[0]][0].values())[0][1][1]
    sample_combined = sample_left + sample_right
    sample_list = [sample_combined, sample_left, sample_right]
    variance_list = []
    try:
        for values in sample_list:
            mean = sum(values)/len(values)
            if len(values) == 1:
                variance = 0
            else:
                sum1 = 0
                for value in values:
                    sum1 += (value - mean)**2
                variance = sum1/(len(values)-1)
            variance_list = variance_list + [variance]
    #try:
        variance_reduction_score = (variance_list[0] - (variance_list[1]*(len(sample_left)/len(sample_combined)) + variance_list[2]*(len(sample_right)/len(sample_combined)))/variance_list[0])

        root_dict_open = root_dict['open']
        dict_of_vars = root_dict_open[tuple(route)]
        print(dict_of_vars)
        print(ssr_[4])
        print(ssr_[5])
        #quit()
        other_attempt = root_dict_open[tuple(route)]
        root_dict_open[tuple(route)] = ssr_[0:2] + [variance_reduction_score]

        dict_of_vars_left = {}
        dict_of_vars_right = {}
        for key in dict_of_vars.keys():
            dict_of_vars_left[key] = [dict_of_vars[key][x] for x in ssr_[4]]
            #print(dict_of_vars_left)
            #print(dict_of_vars_right)
            #print(dict_of_vars)
            dict_of_vars_right[key] = [dict_of_vars[key][x] for x in ssr_[5]]
        left = route + [0]
        right = route + [1]
        if len(left) > 10:
            pass
        root_dict_open[tuple(left)] = dict_of_vars_left
        root_dict_open[tuple(right)] = dict_of_vars_right
        root_dict['open'] = root_dict_open
        complete_MULTI_tree(root_dict=root_dict, route=left)
    except:
        #print(variance_list)
        #print(sample_combined)
        #print(root_dict_open[tuple(route)])
        #print(route)
        #print(root_dict)
        #print(other_attempt)

        values_ = alternate_use[dep_name]
        mean_ = sum(values_)/len(values_)
        root_dict_closed = root_dict['closed']
        root_dict_closed[tuple(route_)] = mean_
        del root_dict_open[tuple(route_)]
        i = 1
        root_dict['open'] = root_dict_open
        root_dict['closed'] = root_dict_closed
        while True:
            #print(route[0:-i + 1])
            try:
                if route_[-i] == 0: # if it is 0, the next option is to change the the same indexed item in the route list to 1
                    if i == 1:
                        route_[-i] = 1
                        break
                    else:
                        #print(route)
                        route_ = route_[0:-i+1]
                        #print(route)
                        route_[-1] = 1
                        break
                else:
                    i += 1
            except IndexError:
                return root_dict
        complete_MULTI_tree(root_dict=root_dict, route=route_)


    

def no_trees(data, n=100):
    dict_of_trees = {}
    for i in range(0, n):
        bootstrapped = bootstrap(data=data)
        #print(bootstrapped)
        #quit()

        root_choice = best_ind_var(bootstrapped=bootstrapped, dep_name='biomarker')
        #print(root_choice)
        #quit()


        multi_choice = initial_multi_root_dict(root_choice=root_choice, bootstrapped=bootstrapped)

        ## SO FAR got the start of the tree, left and right branch, in similar format as to how i made it before with a single variable
        # next steps
        # use best_ind_var() on the branches, then repeat branch splitting with a modified complete_tree()


        #
        #root choice guves you the variable name, and the ssr, and the indexes remixed... 
        # need to put into complete multi tree, with slight diff iterations
        final_tree = complete_MULTI_tree(multi_choice)

        dict_of_trees[n] = final_tree
    return dict_of_trees