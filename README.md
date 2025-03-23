# RandomForest_Project
This is my project writing a Random-Forest based ML algorithm with minimal/no prewritten code. 

Here is an updated (23032025) version of the testing scripts so far. 
  - Currently it is only for continuous variables and in the steps to determine Important Features for the dependent variable you are testing.
    In my case, I am assessing just one of the variables I have generated.
  - Errors: a few, ranging from an index error if one of the branches determines that all of the samples grouped provides the lowest SSRs
    another is that it works once as a TREE, but if you loop it (to generate a FOREST of trees) it appears to 'skip' analysing some of the branches of the trees. No idea why, the code being run is no different to if ran once.
