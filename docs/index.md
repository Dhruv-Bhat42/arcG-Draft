# arcG: Associators for Representation Categories of finite Groups

arcG is a software tool that can compute the F symbols for \(\operatorname{Rep}(G) \) - the representation categories of finite groups. 

For end-user instructions on installation and the use of arcG, please click on the User Guide tab.

arcG consists of a set of GAP scripts:

1. ```G_fusion_finder.g```
2. ``` G_irrep_finder.g```
3. ```label_methods.g```
4. ```ModifiedGS.g```
5. ```main.g```

For a detailed documentation on each file, please click on the corresponding tab. 

# Software Flow Summary
## Input Phase: 
arcG asks the user:

1. to input the name of a finite group (with sample format)
2. to input a range for random integer matrices used in the main algorithm to compute full rank intertwiners between pairs of irreps and their “block diagonal” form 

## Pre-Processing:
Loads character table and irreps, and unitarizing the irreps*. 
*GAP doesn’t have the preset capability to unitarize irreps for every group (we have an exhaustive list of every group of order < 32 for which this fails)
Records and stores runtimes

## Processing
Computes fusion rules via character table
Main Step- Computes full rank intertwiners via “main algorithm” (using invertible random matrices)
Orthogonalizes columns of the matrices of full rank intertwiners via Modified Gram Schmidt 
Extracts trivalent vertices and stores their normalization (i.e. creates a matrix-scalar tuple for every trivalent vertex)
Records and stores runtimes

## Computing F symbols
Creates new data structures to aid in the computation of F symbols (a dictionary of fusion rules, trivalent vertices, triplet fusion outcomes, etc.)
Performs appropriate “diagram stacking” (i.e. tensoring and composing the correct matrices) and extracts the F symbol in 2 parts-
The trace of the final matrix resulting from diagram stacking
The product of all normalizations (i.e. normalization of the F symbol)
Records and stores runtimes

## Output

The output of arcG consists of the following files:
1. An Overview File containing the name of the group G, the GAP ID of G, total number of irreps, whether multiplicities are present, irrep labels, fusion rules, and the number of non-trivial F matrices & F symbols computed. A list of compute times for Preprocessing, Computing Full Rank Intertwiners, Orthogonalization, F Symbols, and Total Compute Time. 

2. List of Unitarized Irreducible Representations of G

3. Main Output File containing all F symbol data of G (with informative headers). Here, arcG sequentially prints F matrices in ascending order (of irrep degree in each index of the F matrix from left to right). It prints each F symbol datum in ascending order of irrep degree in two parts:
- F symbol numerator (computed by diagram stacking)
- F symbol denominator (computed by multiplying normalizations of every trivalent vertex involved in diagram stacking)

[INSERT PICTURES]


