[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# [Decomposable Formulation of Transmission Constraints for Decentralized Power Systems Optimization](https://doi.org/10.1287/ijoc.2022.0326)

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [BSD License](LICENSE).

This repository contains supporting material for the paper 
 
    "Decomposable Formulation of Transmission Constraints for Decentralized Power Systems Optimization" by Álinson S. Xavier, Feng Qiu and Santanu S. Dey.

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2022.0326

https://doi.org/10.1287/ijoc.2022.0326.cd

Below is the BibTex for citing this snapshot of the respoitory.

```
@article{Xavier2024,
  author =        {\'Alinson S. Xavier and Feng Qiu and Santanu S. Dey},
  publisher =     {INFORMS Journal on Computing},
  title =         {Decomposable Formulation of Transmission Constraints for Decentralized Power Systems Optimization},
  year =          {2024},
  doi =           {10.1287/ijoc.2022.0326},
  url =           {https://github.com/INFORMSJoC/2022.0326},
}
```


## Description

One of the most complicating factors in decentralized solution methods for a broad range of power system optimization problems is the modeling of power flow equations. Existing formulations for direct current (DC) power flows either have limited scalability or are very dense and unstructured, making them unsuitable for large-scale decentralized studies. In this work, we present a novel sparsified variant of the injection shift factors formulation, which has a decomposable block-diagonal structure and scales well for large systems. We also propose a decentralized solution method, based on alternating direction multiplier method (ADMM), that efficiently handle transmission line outages in N-1 security requirements. Benchmarks on Multi-Zonal Security-Constrained Unit Commitment problems show that the proposed formulation and algorithm can reliably and efficiently solve interconnection-level test systems, with up to 6,515 buses, with no convergence or numerical issues.


## Requirements

The code was tested on Ubuntu Linux 20.04.3 LTS, with the following software:

- Julia 1.6.7
- IBM® ILOG® CPLEX® Optimization Studio 12.9
- MPICH 3.3.2

## Replicating

1. Install Julia, CPLEX and MPICH

2. Navigate to the project folder and initialize the project:

    ```
    cd Decomposition/
    julia --project=.
    ]dev deps/ADMM deps/unitcommitment.jl
    ]build
    ```

3. Navigate to scripts and run the experiments:

    ```
    cd Decomposition/scripts/
    ./run.sh all
    ```

4. The raw results will be stored in `Decomposition/scripts/logs`

## License

This software is released under the BSD license. See file 'LICENSE' for more information.

