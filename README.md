# pairednetworkbootstrap
R code for implementing paired network inference for test of equality and test of scaling for a wide variety of random graph models. Includes methods proposed in the paper "A Bootstrap-based Inference Framework for Testing Similarity of Paired Networks" by Somnath Bhadra, Kaustav Chakraborty, Srijan Sengupta, and Soumendra Lahiri. Please see the arXiv version of the paper at https://arxiv.org/abs/1911.06869 to use the codes. Codes also include tests proposed in two previous papers, as described below.

Use all_models.R as your source to load the R functions for paired network inference. Description for the functions is given below. Please cite the paper "A Bootstrap-based Inference Framework for Testing Similarity of Paired Networks" by Somnath Bhadra, Kaustav Chakraborty, Srijan Sengupta, and Soumendra Lahiri (arXiv version at https://arxiv.org/abs/1911.06869) if you use the codes.

##############################################################################
T.ase.equal

Description

Computes the p-value of the test for the equality case in rdpg model, proposed by Tang et. al. (2017).

Usage

T.ase.equal(A1, A2, d, B)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
d	Number of parameters for Adjacency Spectral Embedding (ASE)
B	Number of Bootstrap resamples

Output

A array of p-values, one for each egonet. 

##############################################################################

T.ase.scale

Description

Computes the p-value of the test for the scaled case in rdpg model, proposed by Tang et. al. (2017).

Usage

T.ase.scale(A1, A2, d, B)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
d	Number of parameters for Adjacency Spectral Embedding (ASE)
B	Number of Bootstrap resamples

Output

A array of p-values, one for each egonet. 

##############################################################################

T.frob.rdpg

Description

Computes the p-value of the test for the equality case in rdpg model, proposed in this paper (Bhadra et al, 2019+).

Usage

T.frob.rdpg(A1, A2, d, B)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
d	Number of parameters for Adjacency Spectral Embedding (ASE)
B	Number of Bootstrap resamples

##############################################################################

T.scale.rdpg

Description

Computes the p-value of the test for the scaled case in rdpg model, proposed in this paper (Bhadra et al, 2019+).

Usage

T.scale.rdpg(A1, A2, d, B)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
d	Number of parameters for Adjacency Spectral Embedding (ASE)
B	Number of Bootstrap resamples

##############################################################################

T.frob.chunglu

Description

Computes the p-value of the test for the equality case in chung-lu model, proposed in the paper (Bhadra et al, 2019+).

Usage

T.frob.chunglu(A1, A2, B)

Arguements

A1	The first adjacency matrix
A2	The second adjacency matrix
B	Number of Bootstrap resamples

##############################################################################

T.scale.chunglu

Description

Computes the p-value of the test for the scaled case in chung-lu model, proposed in the paper (Bhadra et al, 2019+).


Usage

T.scale.chunglu(A1, A2, B)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
B	Number of Bootstrap resamples

##############################################################################

T.frob.dcbm

Description

Computes the p-value of the test for the equality case in dcbm, proposed in this paper (Bhadra et al, 2019+).

Usage

T.frob.dcbm(A1, A2, k, B)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
k	Number of communities or groups for the dcbm
B	Number of Bootstrap resamples

##############################################################################

T.scale.dcbm

Description

Computes the p-value of the test for the scaled case in dcbm, proposed in this paper (Bhadra et al, 2019+).

Usage

T.scale.dcbm(A1, A2, k, B)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
k	Number of communities or groups for the dcbm
B	Number of Bootstrap resamples

##############################################################################

TWtest

Description

Computes the p-value of the test for the equality case in rdpg model, proposed by Ghoshdastidar and von Luxburg (2018), using cutoff from 
Tracy Widom distribution.

Usage

TWtest(A1, A2, r)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
r	Precomputed communities

##############################################################################

T.frob.pabm

Description

Computes the p-value of the test for the equality case in pabm model, proposed in this paper (Bhadra et al, 2019+).

Usage

T.frob.pabm(A1, A2, k, B)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
k	Number of communities or groups for the pabm
B	Number of Bootstrap resamples

##############################################################################

T.scale.pabm

Description

Computes the p-value of the test for the scaled case in pabm model, proposed in this paper (Bhadra et al, 2019+).

Usage

T.scale.pabm(A1, A2, k, B)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
k	Number of communities or groups for the pabm
B	Number of Bootstrap resamples

##############################################################################

T.frob.pabm.par

Description

Computes the p-value of the test for the equality case in pabm model, proposed in this paper (Bhadra et al, 2019+), using parallel computing to decrease the 
ruuning time of the code by several times.

Usage

T.frob.pabm.par(A1, A2, k, B, cores)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
k	Number of communities or groups for the pabm
B	Number of Bootstrap resamples
cores	Number of cores of the machine one will use to run the code while using parallel computing

##############################################################################

T.scale.pabm.par

Description

Computes the p-value of the test for the scaled case in pabm model, proposed in this paper (Bhadra et al, 2019+), using parallel computing to decrease the 
ruuning time of the code by several times.

Usage

T.scale.pabm.par(A1, A2, k, B, cores)

Arguments

A1	The first adjacency matrix
A2	The second adjacency matrix
k	Number of communities or groups for the pabm
B	Number of Bootstrap resamples
cores	Number of cores of the machine one will use to run the code while using parallel computing

################################################################################
