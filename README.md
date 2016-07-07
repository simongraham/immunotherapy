# Immunotherapy

The code included looks to model the immunotherapy process.

The model created is a cellular automaton, where black cells represent blank space, red cells represent cancer cells, blue cells represent T-cells and green cells represent healthy cells. 

Every iteration, a random number between 0 and 1 is generated. If this number is less than the cancer growth rate, then a randomly chosen cancer cell will grow. Otherwise, all T-cells will move. 

If a T-cell has a cancer cell as its neighbour, there will be a probability that the T-cell will kill the cancer cell. This probability is modelled by an exponential function, where a T-cell is more likely to react with a similar cancer antigen. 

Unlike cancer cells, healthy cells do not possess specific antigens and will be killed by a T-cell with a given small probability. The healthy cell grow rate is density dependent. 

We also model the adoptive immunotherapy process. At a given time, a sample of T-cells are taken, and then are used to repopulate the grid. 

The questions we pose are:

- How does the diversity T-cell clonotypes and cancer cells affect treatment.
- How does the level of cross reactivity affect cancer treatment. Is there a trade off due to the presence of healthy cells?

NB. Cross reactivity denotes the ability for a T-cell to interact with antigens that are not necessarily similar.


