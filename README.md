# Lifted Collocation Integrators for Direct Optimal Control in ACADO Toolkit

This repository contains supporting material for the paper: "Lifted Collocation Integrators for Direct Optimal Control in ACADO Toolkit".
The code allows one to reconstruct the presented numerical results for the optimal control of a chain of masses. The numerical results have been obtained using the open-source ACADO code generation tool on a standard computer, equipped with Intel i7-3720QM processor, and using a 64-bit version of Ubuntu 14.04 and the g++ compiler version 4.8.4.


## General Instructions

First install the Matlab interface of the ACADO Toolkit, which is included as a zip archive '*acado-master.zip*' in this repository. The installation instructions can be found below.

Running the Matlab scripts in the folder chain_mass/ACADO_SQP will generate and compile the C-code for the OCP solver and runs the SQP method based on the corresponding lifted collocation integrator. Note that all these Matlab scripts define a flag '*EXPORT*' and '*COMPILE*' which you can set to zero in case you respectively do not need to export or recompile the solver again (in order to speedup the process of running the simulations).

- By running such a Matlab script, e.g. '*chain_mass/ACADO_SQP/Gauss_Newton/chain_mass.m*', the results are written to a .mat file which can later be used to generate the figures and/or tables from the paper. There are 4 variants of the algorithms which can each be tested both within the Gauss-Newton and the exact Hessian based SQP method:

  - without lifting (MS): '*chain_mass/ACADO_SQP/Gauss_Newton/chain_mass.m*' and '*chain_mass/ACADO_SQP/Exact_Hessian/chain_mass.m*'

  - exact lifting (LC-EN): '*chain_mass/ACADO_SQP/Gauss_Newton/chain_mass_EN.m*' and '*chain_mass/ACADO_SQP/Exact_Hessian/chain_mass_EN.m*'

  - IN lifting (LC-IN): '*chain_mass/ACADO_SQP/Gauss_Newton/chain_mass_IN.m*' and '*chain_mass/ACADO_SQP/Exact_Hessian/chain_mass_IN.m*'

  - INIS lifting (LC-INIS): '*chain_mass/ACADO_SQP/Gauss_Newton/chain_mass_INIS.m*' and '*chain_mass/ACADO_SQP/Exact_Hessian/chain_mass_INIS.m*'
				
- To produce the convergence plots based on the results saved in the .mat files, you can run the Matlab script '*chain_mass/makeFigures.m*'. 

- To produce the tables detailing the computation times based on the saved results, you can run the Matlab script '*chain_mass/ACADO_SQP/makeTables.m*'.


Note that the above mentioned simulation scripts contain multiple parameters related to the optimal control problem formulation and its discretization:

   - Ns: number of integration steps per shooting interval
   - d: number of Gauss-Legendre collocation nodes (d = 4 --> 4-stage Gauss method of order 8, d = 3 --> 3-stage Gauss method of order 6)
   - Nm: number of masses in the chain mass problem, such that the number of states corresponds to nx = (Nm-1)*6  (+1 for the time variable in case of the time optimal formulation)
   - T: the control horizon length as defined in the paper (= 5.0 s for the minimum effort problem, = 1.0 s for the time optimal formulation)
   - N: number of shooting intervals in the multiple shooting discretization
   - wall_pos: this is the position of the wall next to the equilibrium position of the chain mass system, this constraint is imposed on all the free masses as follows:

		for i = 1:Nm-1
			ocp.subjectTo( wall_pos <= ps{i}(2) );
		end

     note that other constraints are included such as the control bounds over the horizon:

		ocp.subjectTo( -10.0 <= u <= 10.0 );

     or the terminal equality constraint for the end state:

		ocp.subjectTo( 'AT_END', states == xN_term );

Each of these scripts exports C-code by calling '*mpc.exportCode(...)*' for an ACADO generated SQP method based on a specific variant of lifted collocation as mentioned above.
This solver is automatically compiled into a MEX file afterwards which can be run with a specific input and output structure of which the most important fields are:

   - input.x0: the initial value constraint in the OCP formulation
   - input.x: the current state trajectory over the horizon
   - input.u: the current control trajectory over the horizon
   - input.y: reference trajectory over the horizon in case of a least squares tracking objective
   - input.yN: the terminal value in the reference trajectory
   - input.W: cost matrix over the horizon in case of a least squares tracking objective
   - input.WN: the terminal cost matrix for the tracking objective

   - output.info: all information which could be interesting to the user regarding the performed SQP step, such as the status flag from the QP solver, the corresponding computation times, etc.
   - output.x: the updated state trajectory over the horizon
   - output.u: the updated control trajectory over the horizon


## ACADO Toolkit: implementation of lifted collocation integrators

All proposed variants of the lifted collocation integrators are implemented in the code generation tool of the ACADO Toolkit and are therefore made available as open-source software.
In this specific implementation, the collocation methods are based on either Gauss-Legendre or Radau IIA points and the proposed Jacobian approximations are based on either Simplified 
or Single Newton-type iterations. The software can be downloaded freely from https://github.com/acado/acado and can be discussed on an active forum https://sourceforge.net/p/acado/discussion/general/.
The ACADO code generation tool is a specific part of this toolkit, which can be used to obtain real-time feasible codes for dynamic optimization on embedded control hardware. 
In particular, it pursues the export of highly efficient C-code based on the RTI scheme for Nonlinear MPC (NMPC). A user friendly MATLAB interface is available which allows one to export, 
compile and use auto generated code in an intuitive way and without direct interaction with C/C++ programming. 

The following integrator classes specifically implement the algorithms presented in this paper, which can be found in the folder '*acado/code_generation/integrators/*':

   - without lifting (MS): 
	the class *ForwardIRKExport* (first order sensitivity propagation) and *SymmetricIRKExport* (including second order sensitivity propagation), 
	respectively defined in the files '*irk_forward_export.cpp*' and '*irk_symmetric_export.cpp*'

   - exact lifted collocation (LC-EN), as well as the inexact INIS lifting (LC-INIS): 
	the class *ForwardLiftedIRKExport* (first order sensitivity propagation) and *SymmetricLiftedIRKExport* (including second order sensitivity propagation), 
	respectively defined in the files '*irk_lifted_forward_export.cpp*' and '*irk_lifted_symmetric_export.cpp*'

   - inexact IN lifting (LC-IN): 
	the class *AdjointLiftedIRKExport* (first and second order sensitivity propagation), in the file '*irk_lifted_adjoint_export.cpp*'


## ACADO Toolkit: installation

(Note that these detailed instructions to install ACADO Toolkit from its Matlab interface can also be found on: http://acado.github.io/matlab_overview.html In addition, there exists an active discussion forum which can be found here: https://sourceforge.net/p/acado/discussion/general/)


To install and use the MATLAB interface you need to have a recent MATLAB version and a C++ compiler installed.
In order to know which compilers are compatible with your Matlab version and on your platform, you can have a look here: http://de.mathworks.com/support/compilers/R2016a/
and to change the default compiler for Matlab, you can follow the instructions here: http://de.mathworks.com/help/matlab/matlab_external/changing-default-compiler.html

- You can then start to compile the source code of ACADO Toolkit from Matlab. For this purpose, you should unzip the archive and go to the main folder of acado_master from Matlab where you then navigate further to the Matlab installation directory:

		cd interfaces/matlab/;

- You are now ready to compile the ACADO interface. This compilation will take several minutes but needs to be done only once. Run “make” in your command window:

		make

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

