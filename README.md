# Lifted Collocation Integrators for Direct Optimal Control in ACADO Toolkit

This repository contains supporting material for the paper: "Lifted Collocation Integrators for Direct Optimal Control in ACADO Toolkit".
The code allows one to reconstruct the presented numerical results for the optimal control of a chain of masses. The numerical results have been obtained using the open-source ACADO code generation tool on a standard computer, equipped with Intel i7-3720QM processor, and using a 64-bit version of Ubuntu 14.04 and the g++ compiler version 4.8.4.

## General Instructions

First install the Matlab interface of the ACADO Toolkit, which is included as a zip archive '*acado-master.zip*' in this repository. The installation instructions can be found below.

Running the Matlab scripts in the folder chain_mass/ACADO_SQP will generate the C-code for the OCP solver and runs the SQP method based on the corresponding lifted collocation integrator. 

- By running such a Matlab script, e.g. '*chain_mass/ACADO_SQP/Gauss_Newton/chain_mass.m*', the results are written to a .mat file which can later be used to generate the figures and/or tables from the paper. There are 4 variants of the algorithms which can each be tested both within the Gauss-Newton and the exact Hessian based SQP method:

  - without lifting (MS): '*chain_mass/ACADO_SQP/Gauss_Newton/chain_mass.m*' and '*chain_mass/ACADO_SQP/Exact_Hessian/chain_mass.m*'

  - exact lifting (LC-EN): '*chain_mass/ACADO_SQP/Gauss_Newton/chain_mass_EN.m*' and '*chain_mass/ACADO_SQP/Exact_Hessian/chain_mass_EN.m*'

  - IN lifting (LC-IN): '*chain_mass/ACADO_SQP/Gauss_Newton/chain_mass_IN.m*' and '*chain_mass/ACADO_SQP/Exact_Hessian/chain_mass_IN.m*'

  - INIS lifting (LC-INIS): '*chain_mass/ACADO_SQP/Gauss_Newton/chain_mass_INIS.m*' and '*chain_mass/ACADO_SQP/Exact_Hessian/chain_mass_INIS.m*'
				
- To produce the convergence plots, you can run the Matlab script '*chain_mass/makeFigures.m*'. 

- To produce the tables detailing the computation times, you can run the Matlab script '*chain_mass/ACADO_SQP/makeTables.m*'.

Optionally, you can also run the scripts in the folder chain_mass/CasADi_IPOPT to obtain a reference solution to the optimal control problem using Ipopt from CasADi (follow the instructions on https://github.com/casadi/casadi/wiki to download the binary for your platform: https://github.com/casadi/casadi/wiki/InstallationInstructions).

## ACADO Toolkit: installation

(Note that these detailed instructions to install ACADO Toolkit from its Matlab interface can also be found on: http://acado.github.io/matlab_overview.html In addition, there exists an active discussion forum which can be found here: https://sourceforge.net/p/acado/discussion/general/)


To install and use the MATLAB interface you need to have a recent MATLAB version and a C++ compiler installed.
In order to know which compilers are compatible with your Matlab version and on your platform, you can have a look here: http://de.mathworks.com/support/compilers/R2016a/
and to change the default compiler for Matlab, you can follow the instructions here: http://de.mathworks.com/help/matlab/matlab_external/changing-default-compiler.html

- You can then start to compile the source code of ACADO Toolkit from Matlab. For this purpose, you should unzip the archive and go to the main folder of acado_master from Matlab where you then navigate further to the Matlab installation directory:

		cd interfaces/matlab/;

- You are now ready to compile the ACADO interface. This compilation will take several minutes but needs to be done only once. Run “make” in your command window:

		make clean all;

- You will see:

		Making ACADO...

- and after a while when the compilation is finished:

		ACADO successfully compiled.
		Needed to compile XYZ file(s).

		If you need to restart Matlab, run this make file again
		to set all paths or run savepath in your console to
		save the current search path for future sessions.

- ACADO has now been compiled. As the text indicated every time you restart MATLAB you need to run “make” again to set all paths. When running “make” again no new files need to be compiled and the process will only take a few seconds. However, it is easier to save your paths for future Matlab session. Do so by running “savepath” in your command window (this step is optional).

		savepath;
