#!/usr/bin/python3


# Install Packages
from numpy.random import seed
from numpy.random import normal
from matplotlib.pyplot import hist
from matplotlib.pyplot import scatter
import matplotlib.pyplot as plt
from one_variable_tree_functions import *
from random import randint
from multi_variable_tree_functions import *


# MORE variables!
sigmoid, random_, effective = more_variables(30) # this is just for the testing
#print(sigmoid) # bioavailability
#print(random_) # biomarker that you check (in this case dose changes nothing)
#print(effective) # patient still has disease or not (1 is disease free)

independent_variable, dependent_variable = stable_normally_distributed_plot(1000, 30)
#print(independent_variable) # mg/kg dose

# TIME TO BOOTSTRAP - my code will assume as a dictionary with keys as column titles
data = {'bioavailability': sigmoid,
'other_biomarker': random_,
'treated': effective,
'mgkg': independent_variable,
'biomarker': dependent_variable}
# theoretically this is what the data should look like when bootstrapped

forest = no_trees(data=data, n=1)


print(forest)









# DONE Now need to turn it into a random forest
# maybe do the feature importance first
