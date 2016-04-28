%
%    This file was auto-generated by ACADO Code Generation Tool.
%    
%    ACADO Code Generation tool is a sub-package of ACADO toolkit --
%    A Toolkit for Automatic Control and Dynamic Optimization.
%    Copyright (C) 2008-2014 by Boris Houska, Hans Joachim Ferreau,
%    Milan Vukov, Rien Quirynen, KU Leuven.
%    Developed within the Optimization in Engineering Center (OPTEC)
%    under supervision of Moritz Diehl. All rights reserved.
%    
%    ACADO Toolkit is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%    
%    ACADO Toolkit is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%    
%    You should have received a copy of the GNU Lesser General Public
%    License along with ACADO Toolkit; if not, write to the Free Software
%    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%    


function make_rien_solver( name )

	% Output file name, and also function name
	if (nargin > 0)
		fileOut = name;
	else
		fileOut = 'acado_solver';
	end;
		
	% Root folder of code generation
	CGRoot = '.';	
	
	% qpOASES embedded source files
	qpOASESSources = [ ...
		'CGRoot/qpoases/SRC/Bounds.cpp ' ...
		'CGRoot/qpoases/SRC/Constraints.cpp ' ...
		'CGRoot/qpoases/SRC/CyclingManager.cpp ' ...
		'CGRoot/qpoases/SRC/Indexlist.cpp ' ...
		'CGRoot/qpoases/SRC/MessageHandling.cpp ' ...
		'CGRoot/qpoases/SRC/QProblem.cpp ' ...
		'CGRoot/qpoases/SRC/QProblemB.cpp ' ...
		'CGRoot/qpoases/SRC/SubjectTo.cpp ' ...
		'CGRoot/qpoases/SRC/Utils.cpp ' ...
		'CGRoot/qpoases/SRC/EXTRAS/SolutionAnalysis.cpp ' ...
		];
		
	% Auto-generated files
	CGSources = [ ...
		'CGRoot/rien_solver_mex.c ' ...
		'CGRoot/acado_solver.c ' ...
		'CGRoot/acado_integrator.c ' ...
		'CGRoot/acado_auxiliary_functions.c ' ...
		'CGRoot/acado_qpoases_interface.cpp ' ...
		];
	if (nargin > 1)
		CGSources = [CGSources extern];
	end
		
	% Adding additional linker flags for Linux
	ldFlags = '';
	if (isunix() && ~ismac())
		ldFlags = '-lrt';
    elseif (ispc)
        ldFlags = '-DWIN32';
	end;

	% Recipe for compilation
	CGRecipe = [ ...
		'mex -O' ...
		' -I. -I''CGRoot'' -I''CGRoot/qpoases'' -I''CGRoot/qpoases/INCLUDE'' -I''CGRoot/qpoases/SRC''' ...
		' ldFlags' ...
		' -D__MATLAB__ -DQPOASES_CUSTOM_INTERFACE="acado_qpoases_interface.hpp" -O CGSources qpOASESSources -output %s.%s' ...
		];

% Compilation
qpOASESSources = regexprep(qpOASESSources, 'CGRoot', CGRoot);
CGSources = regexprep(CGSources, 'CGRoot', CGRoot);

CGRecipe = regexprep(CGRecipe, 'CGRoot', CGRoot);
CGRecipe = regexprep(CGRecipe, 'CGSources', CGSources);
CGRecipe = regexprep(CGRecipe, 'qpOASESSources', qpOASESSources);
CGRecipe = regexprep(CGRecipe, 'ldFlags', ldFlags);

% disp( sprintf( CGRecipe, fileOut, mexext ) ); 
fprintf( 'compiling... ' );
eval( sprintf(CGRecipe, fileOut, mexext) );
fprintf( ['done! --> ' fileOut '.' mexext '\n'] );
