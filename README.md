# R-tree directional queries
This is a CMake project that allows one to perform directional (and linear) top-k queries using an R-Tree as index.

## Dependencies
This project depends on the [NLopt](https://nlopt.readthedocs.io/en/latest/) library for solving non-linear problems required to execute the directional queries.

## Features

- Perform sequential linear and directional top-k queries given a CSV file containing the dataset.
- Perform linear and directional top-k queries using an R-Tree index given a CSV file containing the dataset, a file containing the values of k to take into account and a file contaning the set of query weights of the queries to perform.
- Perform directional queries using the *TIGHT* method (I/O optimality of accesses).
- Perform directional queries using the *LOOSE* method (no non-linear problems to solve at the cost of higher I/O accesses).
- Track the performances of the queries (execution time, memory accesses, problems solved...).

## How to use

In the main function, edit the following variables:
-  `datasetPath` : path of the CSV file containing the dataset to process.
-  `kValuesFile` : path of the txt file containing the set of k values for the queries to perform.
-  `queriesFile` : path of the txt file containing the set of query weights for the queries to perform.

Edit the following macros:
-  `DIM` : number of dimensions of the dataset to consider.
-  `BETA` : value of beta to use for the directional query.
-  `TREE METHOD` : set it to *TIGHT* or *LOOSE* .

## Sample files
For the sake of simplicity, we provide some examples of files used and the formatting required in order to properly run the project:

- `datasetPath` : CSV file containing the dataset to process. NOTE: values must be normalized in range [0,1] - [Sample dataset](https://github.com/stefuma19/RtreeV2/tree/master/datasets/cor_neg_1k_2.csv).
-  `kValuesFile` : TXT file containing the set of k values for the queries to perform - [Sample k](https://github.com/stefuma19/RtreeV2/blob/master/utilities/k.txt).
-  `queriesFile` : TXT file containing the set of query weights for the queries to perform - [Sample queries](https://github.com/stefuma19/RtreeV2/blob/master/queries/2d.txt).

## Output
The output consists of a CSV file that will be located in the project root folder. This file, `stats.csv`,  contains statistics about the dataset used, the corresponding R-Tree and the mean of the results of the (linear and directional) queries performed.

## Output explanation
This section presents the explanation of the columns of the output file.
- `n` : Cardinality of the dataset.
- `d` : Dimensionality of the dataset.
- `k` : Top-*k* value considered.
- `Rtree boxes` : Number of boxes (or nodes) of the built Rtree.
- `Rtree leaves` : Number of leaf nodes of the built Rtree.
- `Method` : Method used to perform the directional queries (Either *TIGHT* or *LOOSE*).
- `Linear RTree time [us]` : Time taken, on average, to perform the linear top-k queries with the value of k and the queries specified.
- `Linear Rtree numPoint` : Number of points accessed, on average, over the linear top-k queries performed.
- `Linear Rtree numLeaves` : Number of leaf nodes accessed of the R-tree, on average, over the linear top-k queries performed.
- `Linear Rtree numBoxes` : Number of nodes (or boxes) accessed of the R-tree, on average, over the linear top-k queries performed.
- `Linear Rtree Bounds Computed` : Number of bounds computed, on average, over the linear top-k queries performed to decide if we can stop visiting a subtree.
- `Linear Rtree Bounds Computation Time [us]` : Time taken to solve the previously described problems, on average, over the linear top-k queries performed.
- `Linear Rtree Scores` : Number of scores computed (it should coincide with `Linear Rtree numPoint`).
- `Linear Rtree Scores Computation Time [us]` : Time taken to solve the previously described scores, on average, over the linear top-k queries performed.
- `Directional RTree time [us]` : Time taken, on average, to perform the directional top-k queries with the value of k and the queries specified.
- `Directional Rtree numPoint` : Number of points accessed, on average, over the directional top-k queries performed.
- `Directional Rtree numLeaves` : Number of leaf nodes accessed of the R-tree, on average, over the directional top-k queries performed.
- `Directional Rtree numBoxes` : Number of nodes (or boxes) accessed of the R-tree, on average, over the directional top-k queries performed.
- `Directional Rtree Problems` : Number of non-linear problems computed, on average, over the directional top-k queries performed to decide if we can stop visiting a subtree.
- `Directional Rtree ProblemTime [us]` : Time taken to solve the previously described problems, on average, over the directional top-k queries performed.
- `Directional Rtree Scores` : Number of directional query scores computed (it should coincide with `Directional Rtree numPoint`).
- `Directional Rtree Scores Computation Time [us]` : Time taken to solve the previously described scores, on average, over the directional top-k queries performed.
- `Linear RTtree Bounds Computed For Directional` : Number of bounds computed for the directional query with the *LOOSE* method.
- `Linear RTtree Bounds Computed For Directional Computation Time [us]` : Time taken to solve the previously described bounds, on average, over the directional top-k queries performed with the *LOOSE* method. 


