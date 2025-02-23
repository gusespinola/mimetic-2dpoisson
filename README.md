# mimetic-2dpoisson
Two-dimensional discretization of the Poisson equation

# Two-dimensional discretization of the Poisson equation by using high-order mimetic operators

This repository contains the necessary files to reproduce the results from [1].

## Requirements

To reproduce the results, you will need the Mimetic Operator Library Enhanced[(MOLE)](https://github.com/csrc-sdsu/mole) (as forked [here](https://github.com/gusespinola/mole)) and the Krylov Subspace-Based Adaptive Solvers [(KrySBAS)](https://github.com/nidtec-una/krysbas-dev/) library (release [0.3.1](https://github.com/nidtec-una/krysbas-dev/releases/tag/v0.3.1)).

## Reproducing the results

To reproduce the results, you should run the file `minimal_poisson2D.mat`.

If you wish to generate plots for convergence and evolution of the parameter m, and/or skip some computations, you can activate the corresponding control flags embedded in the source code.

## Citing

If you use (either all or a part of) the code present in this repository, I ask you to cite [1].

## References

[1] **Espínola, G. E., Cabral, J. C., Schaerer, C. E., & Varela, J.** (2025). *A numerical investigation on iterative methods for the two-dimensional Poisson equation discretized with high-order mimetic operators*. Proceeding Series of the Brazilian Society of Computational and Applied Mathematics, 11(1), 1-2. URL: [https://proceedings.sbmac.org.br/sbmac/article/view/4723](https://proceedings.sbmac.org.br/sbmac/article/view/4723)
