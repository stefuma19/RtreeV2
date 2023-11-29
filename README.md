# R-tree directional queries
This is a CMake project that allows one to perform directional (and linear) top-k queries using an R-Tree as index.

## Dependencies
This project depends on the [NLopt](https://nlopt.readthedocs.io/en/latest/) library for solving non-linear problems required to execute the directional queries.

## Features

- Perform sequential linear and directional top-k queries given a CSV file containing the dataset.
- Perform linear and directional top-k queries using an R-Tree index given a CSV file containing the dataset, a file containing the values of k to take into account and a file contaning the set of query weights of the queries to perform.
- Perform directional queries using the TIGHT method (I/O optimality of accesses).
- Perform directional queries using the LOOSE method (no non-linear problems to solve at the cost of higher I/O accesses).
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