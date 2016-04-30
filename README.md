# liftedCollocation

This repository contains supporting material for the paper: "Lifted Collocation Integrators for Direct Optimal Control in ACADO Toolkit".
The code allows one to reconstruct the presented numerical results for the optimal control of a chain of masses. The numerical results have been obtained using the open-source ACADO code generation tool on a standard computer, equipped with Intel i7-3720QM processor, and using a 64-bit version of Ubuntu 14.04 and the g++ compiler version 4.8.4.

## Instructions

First install the Matlab interface of the ACADO Toolkit which is a submodule of this repository (follow the instructions on http://acado.github.io/matlab_overview.html to install ACADO from its Matlab interface).

Running the Matlab scripts in the folder chain_mass/ACADO_SQP will generate the C-code for the OCP solver and runs the SQP method based on the corresponding lifted collocation integrator.

As a reference, you can also run the scripts in the folder chain_mass/CasADi_IPOPT to solve the optimal control problem using Ipopt from CasADi (follow the instructions on https://github.com/casadi/casadi/wiki to download the binary).
