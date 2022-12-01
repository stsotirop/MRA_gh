# MRA_gh

## Input parameters:
Max Lengths:
- Maximum lengths for the different orders. The length of the list is the number of the orders that the algorithm will run
- List of numbers

Ratio Limit:
- Ratio of the limit for forces evaluation
- Number

Stiffness:
- Springs stiffness for the n-order method
- Number

Mass:
- Mass of each node
- Number

Time Step: 
- Time step for the n-order method
- Number

Max Iterations: 
- Number of maximum iterations for the n-order method
- Number

ConPointS: 
- Start of constrained points
- List of Point3d (Grasshopper object)

ConPointE: 
- End of constrained points
- List of Point3d (Grasshopper object)

UnConPointS: 
- Start of unconstrained points
- List of Point3d (Grasshopper object)

UnConPointE: 
- End of unconstrained points
- List of Point3d (Grasshopper object)

Stiffness RM:
- Springs stiffness for the method of repulsive masses
- Number

Time Step RM: 
- Time step for the method of repulsive masses
- Number

Max Iterations RM: 
- Number of maximum iterations for the method of repulsive masses
- Number

Run MRA RM: 
- Run or not the method of repulsive masses
- Boolean

## Output parameters:
Lasche:
- Lasche ropes for each order. 
- Data tree of List of integers.
Num Loose:
- Number of loose ropes for each order. 
- List of integers.
Mean Loose:
- Average length of loose ropes for each order. 
- List of numbers.
Min Loose:
- Minimum length of loose ropes for each order. 
- List of numbers.
Num over:
- Number of over ropes for each order. 
- List of integers.
Mean over:
- Average length of over ropes for each order. 
- List of numbers.
Max over:
- Maximum length of loose ropes for each order. 
- List of numbers.
Per Dif:
- Percentage (%) of ropes that are different from lmax for each order. 
- List of numbers.
The outputs: “Connectivity”, “Tolerance”, “Final Points”, “Lfinal” are necessary data for the representation of the structure after the end of each order. The user must connect these outputs to the equivalent inputs of the “MRAplot” component.

# MRAplot

## Input parameters:
Max Lengths (same with MRA):
- Maximum lengths for the different orders. The length of the list is the number of the orders that the algorithm will run. 
- List of numbers.
Display Order:
- Defines the structure of the order that the user wants to produce. If the user selects a number bigger than the orders that the method ran, the component is not running.
- Integer.

## Output parameters:
Order Points:
- The coordinate of the points for the selected order. 
- List of Point3d (Grasshopper object)
Ropes:
- The final ropes for the selected order. 
- List of Line (Grasshopper object)
