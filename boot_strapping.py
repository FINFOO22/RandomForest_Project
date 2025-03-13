#!/usr/bin/python3

# Bootstrapping
# lots of replicate samples (same no. as original sample size, same sample can be picked more than once)


# Where does random forest come into play? 
# To reduce over fitting, colinearisation of trees (because some variables may always picked first at top of tree)
# only use RANDOM subset of variables in each tree to create the tree predict used to predict the dependent variable
# THIS IS DONE AT EACH STEP, so a new set of random variables are picked to predict the dependent variable (ignoring the first one that was used at the top of the tree)


