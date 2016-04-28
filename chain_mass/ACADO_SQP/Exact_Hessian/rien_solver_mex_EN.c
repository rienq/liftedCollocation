/*
 *    This file was auto-generated by ACADO Code Generation Tool.
 *    
 *    ACADO Code Generation tool is a sub-package of ACADO toolkit --
 *    A Toolkit for Automatic Control and Dynamic Optimization.
 *    Copyright (C) 2008-2014 by Boris Houska, Hans Joachim Ferreau,
 *    Milan Vukov, Rien Quirynen, KU Leuven.
 *    Developed within the Optimization in Engineering Center (OPTEC)
 *    under supervision of Moritz Diehl. All rights reserved.
 *    
 *    ACADO Toolkit is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *    
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *    
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *    
 */


/*
input
{
	control: optional
		0: init once, and run preparation and feedback; default behaviour
		1: initialize
		2: preparation
		3: feedback
		4: shift
	x
	u
	mu
	od
	x0: depends on the type of an OCP
	xAC: optional
	SAC: optional
	shifting: optional
	{
		strategy:
			1: use xEnd
			2: integrate
		xEnd
		uEnd
	}
	initialization: optional
		1: initialize by a forward simulation
		else: do nothing
}

output
{
	x
	u
	mu
	xAC: optional
	SAC: optional
	info
	{
		status
		cpuTime
		kktValue
		objValue
		nIterations: works only for qpOASES
	}	
}
*/

/** MEX interface for the ACADO OCP solver
 *
 *  \author Milan Vukov, milan.vukov@esat.kuleuven.be
 *
 *  Credits: Alexander Domahidi (ETHZ), Janick Frasch (KUL, OVGU)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "mex.h"
#include "acado_common.h"
#include "acado_auxiliary_functions.h"

#define FREE( mem ) { if( mem ) { mxFree( mem ); mem = NULL; } }

/* Define number of outputs */
#define NOO 4

#if ACADO_NXA > 0
#define NOO_2 NOO + 1
#else
#define NOO_2 NOO
#endif

#if ACADO_USE_ARRIVAL_COST == 1
#define NOO_3 NOO_2 + 2
#else
#define NOO_3 NOO_2
#endif

#if ACADO_COMPUTE_COVARIANCE_MATRIX == 1
#define NOO_4 NOO_3 + 1
#else
#define NOO_4 NOO_3
#endif

/** Instance of the user data structure. */
ACADOvariables acadoVariables;
/** Instance of the private workspace structure. */
ACADOworkspace acadoWorkspace;


double simTime, qpTime, condenseTime, regularizeTime;

extern int acado_modelSimulation(  );
extern int acado_evaluateObjective(  );
extern int acado_regularizeHessian(  );
extern int acado_condensePrep(  );

extern int acado_condenseFdb(  );
extern int acado_solve(  );
extern int acado_expand(  );

int rien_preparationStep(  )
{
int ret;
acado_timer tmr;

acado_tic( &tmr );
ret = acado_modelSimulation();
simTime = acado_toc( &tmr );

acado_evaluateObjective(  );

acado_tic( &tmr );
acado_regularizeHessian(  );
regularizeTime = acado_toc( &tmr );

acado_tic( &tmr );
acado_condensePrep(  );
condenseTime = acado_toc( &tmr );
return ret;
}

int rien_feedbackStep(  )
{
int tmp;
acado_timer tmr;

acado_tic( &tmr );
acado_condenseFdb(  );
condenseTime += acado_toc( &tmr );

/* printMatrix( "hess", acadoWorkspace.H, 73, 73 );
   printMatrix( "grad", acadoWorkspace.g, 73, 1 ); */

acado_tic( &tmr );
tmp = acado_solve( );
qpTime = acado_toc( &tmr );

acado_tic( &tmr );
acado_expand(  );
condenseTime += acado_toc( &tmr );
return tmp;
}


/** A bit more advanced printing function. */
void mexErrMsgTxtAdv(	char* string,
						...
						)
{
	static char buffer[ 128 ];
	
	va_list printArgs;
	va_start(printArgs, string);
	
	vsprintf(buffer, string, printArgs);
	va_end( printArgs );

	mexErrMsgTxt( buffer );
}

/** A simple helper function. */
void printMatrix(	const char* name,
					real_t* mat,
					unsigned nRows,
					unsigned nCols
					)
{
    unsigned r, c;
    mexPrintf("%s: \n", name);
    for (r = 0; r < nRows; ++r)
    {
        for(c = 0; c < nCols; ++c)
            mexPrintf("\t%e", mat[r * nCols + c]);
        mexPrintf("\n");
    }
}

/** A function for copying data from MATLAB to C array. */
int getArray(	const unsigned mandatory,
				const mxArray* source,
				const int index,
				const char* name,
				real_t* destination,
				const unsigned nRows,
				const unsigned nCols
				)
{
	mxArray* mxPtr = mxGetField(source, index, name);
	unsigned i, j;
	double* dPtr;
	
	if (mxPtr == NULL)
	{
		if ( !mandatory )
			return -1;
		else
			mexErrMsgTxtAdv("Field %s not found.", name);
	}

    if ( !mxIsDouble( mxPtr ) )
		mexErrMsgTxtAdv("Field %s must be an array of doubles.", name);

    if (mxGetM( mxPtr ) != nRows || mxGetN( mxPtr ) != nCols )
		mexErrMsgTxtAdv("Field %s must be of size: %d x %d.", name, nRows, nCols);

	dPtr = mxGetPr( mxPtr );
	
	if (destination == NULL)
		destination = (real_t*)mxCalloc(nRows * nCols, sizeof( real_t ));

	if (nRows == 1 && nCols == 1)
		*destination = *dPtr;
	else
		for (i = 0; i < nRows; ++i)
			for (j = 0; j < nCols; ++j)
				destination[i * nCols + j] = (real_t)dPtr[j * nRows + i];
			
	return 0;
}

void setArray( 	mxArray* destination,
				const int index,
				const char* name,
				real_t* source,
				const unsigned nRows,
				const unsigned nCols
				)
{
	mxArray* mxPtr = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
	double* dPtr = mxGetPr( mxPtr );
	unsigned i, j;
	
	if (nRows == 1 && nCols == 1)
		*dPtr = *source;
	else
		for (i = 0; i < nRows; ++i)
			for(j = 0; j < nCols; ++j)
				dPtr[j * nRows + i] = (double)source[i * nCols + j];

	mxSetField(destination, index, name, mxPtr);
}

/** The MEX interface function. */
void mexFunction(	int nlhs,
					mxArray *plhs[],
					int nrhs,
					const mxArray *prhs[]
					)
{
	static unsigned initialized = 0;
	unsigned ctrl;
	int ctrlIndex, i, j;
	unsigned strategy;
	unsigned initType;
	real_t* xEnd = NULL;
	real_t* uEnd = NULL;
	const mxArray* src = prhs[ 0 ];
	
	const char *infoNames[ 9 ] = {"status", "cpuTime", "simTime", "qpTime", "condenseTime", "regularizeTime", "kktValue", "objValue", "nIterations"};
	mxArray* info;
	double status, cpuTime, kktValue, objValue;
	double tmp[ 1 ];
	mxArray* shPtr;
	acado_timer tmr;
	double nIterations = 0;
	
	const char *outNames[ NOO_4 ];
	outNames[ 0 ] = "info";
	outNames[ 1 ] = "x";
	outNames[ 2 ] = "u";
	outNames[ 3 ] = "mu";
#if ACADO_NXA
	outNames[ NOO ] = "z";
#endif	
		
#if ACADO_USE_ARRIVAL_COST == 1
	outNames[ NOO_2 ] = "xAC";
	outNames[NOO_2  + 1] = "SAC";
#endif
#if ACADO_COMPUTE_COVARIANCE_MATRIX == 1
	outNames[ NOO_3 ] = "sigmaN";
#endif
	
	if (nrhs != 1)
		mexErrMsgTxt(
			"This function requires exactly one input: a structure with parameters.");
			
	if (nlhs != 1)
		mexErrMsgTxt(
			"This function returns one output.");
			
	if( !mxIsStruct( src ) )
		mexErrMsgTxt("The function argument must be a structure.");
	
	/* Get the control flag. */
	if (getArray(0, src, 0, "control", tmp, 1, 1) == 0)
		ctrl = (unsigned)tmp[ 0 ];
	else
		ctrl = 0;
		
	/* Get the initialization flag. */
	if (getArray(0, src, 0, "initialization", tmp, 1, 1) == 0)
		initType = (unsigned)tmp[ 0 ];
	else
		initType = 0;
		
	/* Copy MATLAB arrays to C arrays. */
	getArray(1, src, 0, "x",  acadoVariables.x, ACADO_N + 1, ACADO_NX);
	getArray(1, src, 0, "u",  acadoVariables.u, ACADO_N, ACADO_NU);
	getArray(1, src, 0, "mu", acadoVariables.mu, ACADO_N, ACADO_NX);
	
#if ACADO_NXA	
	getArray(1, src, 0, "z",  acadoVariables.z, ACADO_N, ACADO_NXA);
#endif	
	
#if ACADO_NOD
	getArray(1, src, 0, "od", acadoVariables.od, ACADO_N + 1, ACADO_NOD);
#endif

#if ACADO_INITIAL_STATE_FIXED
	getArray(1, src, 0, "x0", acadoVariables.x0, ACADO_NX, 1);
#endif /* ACADO_INITIAL_STATE_FIXED */

#if (ACADO_HARDCODED_CONSTRAINT_VALUES == 0) && ( (ACADO_QP_SOLVER == ACADO_QPOASES) || (ACADO_QP_SOLVER == ACADO_QPOASES3) )

	if (!initialized)
	{
		acado_initializeSolver();
	}
	
	/* Bounds */
#if ACADO_INITIAL_STATE_FIXED == 1
	getArray(1, src, 0, "lbValues", acadoVariables.lbValues, ACADO_N * ACADO_NU, 1);
	getArray(1, src, 0, "ubValues", acadoVariables.ubValues, ACADO_N * ACADO_NU, 1);
#else
	getArray(1, src, 0, "lbValues", acadoVariables.lbValues, ACADO_NX + ACADO_N * ACADO_NU, 1);
	getArray(1, src, 0, "ubValues", acadoVariables.ubValues, ACADO_NX + ACADO_N * ACADO_NU, 1);
#endif /* ACADO_INITIAL_STATE_FIXED == 0 */

#if QPOASES_NCMAX > 0
	/* Affine constraints */
	getArray(1, src, 0, "lbAValues", acadoVariables.lbAValues, QPOASES_NCMAX, 1);
	getArray(1, src, 0, "ubAValues", acadoVariables.ubAValues, QPOASES_NCMAX, 1);
#endif /* QPOASES_NCMAX > 0 */

#endif /* (ACADO_HARDCODED_CONSTRAINT_VALUES == 0) && ( (ACADO_QP_SOLVER == ACADO_QPOASES) || (ACADO_QP_SOLVER == ACADO_QPOASES3) ) */

#if (ACADO_QP_SOLVER == ACADO_QPDUNES)
	if (!initialized) {
		acado_initializeSolver();
	}
#endif

#if ACADO_USE_ARRIVAL_COST == 1
	getArray(1, src, 0, "xAC", acadoVariables.xAC, ACADO_NX, 1);
	getArray(1, src, 0, "SAC", acadoVariables.SAC, ACADO_NX, ACADO_NX);
    getArray(1, src, 0, "WL", acadoVariables.WL, ACADO_NX, ACADO_NX);
#endif

	/* Shifting strategy */
	shPtr = mxGetField(src, 0, "shifting");
	if (shPtr != NULL)
	{
		if( !mxIsStruct( shPtr ) )
			mexErrMsgTxt("Field \"shifting\" must be defined with a structure.");
		
		/* Get the shifting strategy flag */
		getArray(1, shPtr, 0, "strategy", tmp, 1, 1);
		strategy = (unsigned)tmp[ 0 ];
		
		if (strategy > 2)
			mexErrMsgTxt("Valid options for the shifting strategy are 1 or 2.");
	
		getArray(0, shPtr, 0, "xEnd", xEnd, ACADO_NX, 1);
		getArray(0, shPtr, 0, "uEnd", uEnd, ACADO_NU, 1);
	}
	else
		strategy = 0;
		
	acado_tic( &tmr );
	
	/* Call solver */
	switch ( ctrl )
	{
		case 0:
			/* Simple operational mode. Run one RTI with optional shifting. */
			
			if ( !initialized )
			{
				memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
			
#if ACADO_HARDCODED_CONSTRAINT_VALUES == 1
				acado_initializeSolver();
#endif /* ACADO_HARDCODED_CONSTRAINT_VALUES == 1 */
                
                for( i = 0; i < ACADO_N*ACADO_RK_NIS; i++ ) {
                    for( j = 0; j < ACADO_RK_NSTAGES*(ACADO_NX+ACADO_NXA); j++ ) {
                        acadoWorkspace.rk_A_traj[i*ACADO_RK_NSTAGES*(ACADO_NX+ACADO_NXA)*ACADO_RK_NSTAGES*(ACADO_NX+ACADO_NXA)+j*ACADO_RK_NSTAGES*(ACADO_NX+ACADO_NXA)+j] = 1.0;
                    }
                }
				
				if (initType == 1)
				{
					acado_initializeNodesByForwardSimulation();
				}
				
#if ACADO_USE_ARRIVAL_COST == 1 
                	acado_updateArrivalCost( 1 );
#endif /* ACADO_USE_ARRIVAL_COST == 1 */
				
				initialized = 1;
			}
			else if (strategy == 1 || strategy == 2)
			{
#if ACADO_USE_ARRIVAL_COST == 1 
                acado_updateArrivalCost( 0 );
#endif /* ACADO_USE_ARRIVAL_COST == 1 */
				
				acado_shiftStates(strategy, xEnd, uEnd);
				acado_shiftControls(uEnd);
			}
			
			rien_preparationStep();
			
			status = (double)rien_feedbackStep();
			
			kktValue = acado_getKKT();
			objValue = acado_getObjective();

#if ( (ACADO_QP_SOLVER == ACADO_QPOASES) || (ACADO_QP_SOLVER == ACADO_QPOASES3) )
			nIterations = (double)acado_getNWSR();
#endif /* ( (ACADO_QP_SOLVER == ACADO_QPOASES) || (ACADO_QP_SOLVER == ACADO_QPOASES3) ) */
			
			break;
		
		case 1:
			/* Initialize */
			
			memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
			
			acado_initializeSolver();
				
			if (initType == 1)
			{
				acado_initializeNodesByForwardSimulation();
			}
			
#if ACADO_USE_ARRIVAL_COST == 1 
                acado_updateArrivalCost( 1 );
#endif /* ACADO_USE_ARRIVAL_COST == 1 */			
			
			break;
		
		case 2:
			/* Preparation step */
			
			rien_preparationStep();
			
			break;
		
		case 3:
			/* Feedback step */
			
			status = (double)rien_feedbackStep();
			
			kktValue = acado_getKKT();
			objValue = acado_getObjective();
			
#if ( (ACADO_QP_SOLVER == ACADO_QPOASES) || (ACADO_QP_SOLVER == ACADO_QPOASES3) )
			nIterations = (double)acado_getNWSR();
#endif /* ( (ACADO_QP_SOLVER == ACADO_QPOASES) || (ACADO_QP_SOLVER == ACADO_QPOASES3) ) */
			
			break;
		
		case 4:
			/* Shifting */
			
#if ACADO_USE_ARRIVAL_COST == 1 
                acado_updateArrivalCost( 0 );
#endif /* ACADO_USE_ARRIVAL_COST == 1 */			
			
			acado_shiftStates(strategy, xEnd, uEnd);
			acado_shiftControls( uEnd );
			
			break;
			
		default:
			/* Return an error */
			mexErrMsgTxt("Unknown control code.");
	}
	
	cpuTime = acado_toc( &tmr );
	
	/* Prepare return argument */
	
	plhs[ 0 ] = mxCreateStructMatrix(1, 1, NOO_4, outNames);
		
	setArray(plhs[ 0 ], 0, "x", acadoVariables.x, ACADO_N + 1, ACADO_NX);
	setArray(plhs[ 0 ], 0, "u", acadoVariables.u, ACADO_N, ACADO_NU);
	setArray(plhs[ 0 ], 0, "mu", acadoVariables.mu, ACADO_N, ACADO_NX);
#if ACADO_NXA > 0
	setArray(plhs[ 0 ], 0, "z", acadoVariables.z, ACADO_N, ACADO_NXA);
#endif	
		
#if ACADO_USE_ARRIVAL_COST == 1
	setArray(plhs[ 0 ], 0, "xAC", acadoVariables.xAC, ACADO_NX, 1);
	setArray(plhs[ 0 ], 0, "SAC", acadoVariables.SAC, ACADO_NX, ACADO_NX);
#endif /* ACADO_USE_ARRIVAL_COST */

#if ACADO_COMPUTE_COVARIANCE_MATRIX == 1
	setArray(plhs[ 0 ], 0, "sigmaN", acadoVariables.sigmaN, ACADO_NX, ACADO_NX);
#endif /* ACADO_COMPUTE_COVARIANCE_MATRIX */

	/* Create the info structure. */
	info = mxCreateStructMatrix(1, 1, 9, infoNames);
		
	setArray(info, 0, "status", &status, 1, 1);
	setArray(info, 0, "cpuTime", &cpuTime, 1, 1);
	setArray(info, 0, "simTime", &simTime, 1, 1);
	setArray(info, 0, "qpTime", &qpTime, 1, 1);
	setArray(info, 0, "condenseTime", &condenseTime, 1, 1);
	setArray(info, 0, "regularizeTime", &regularizeTime, 1, 1);
	setArray(info, 0, "kktValue", &kktValue, 1, 1);
	setArray(info, 0, "objValue", &objValue, 1, 1);
	
#if ( (ACADO_QP_SOLVER == ACADO_QPOASES) || (ACADO_QP_SOLVER == ACADO_QPOASES3) )
	setArray(info, 0, "nIterations", &nIterations, 1, 1);
#endif /* ( (ACADO_QP_SOLVER == ACADO_QPOASES) || (ACADO_QP_SOLVER == ACADO_QPOASES3) ) */
		
	mxSetField(plhs[ 0 ], 0, "info", info);
	
	/* Cleanup of the allocated memory */
	FREE( xEnd );
	FREE( uEnd );
}
