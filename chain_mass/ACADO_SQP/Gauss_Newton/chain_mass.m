%% Gauss-Newton based SQP implementation using collocation without lifting
clc;
clear all;
close all;

method = '';

% -------------
% Ns: number of integration steps per shooting interval
Ns = 3; 

% -------------
% d: number of collocation nodes
% d = 4 --> 4-stage Gauss method of order 8
% d = 3 --> 3-stage Gauss method of order 6
d = 4;

% -------------
% Nm: number of masses in the chain
for Nm = 3:8
    
    global ACADO_;
    clearvars -global
    
% Environment
g = 9.81;     % [N/kg]
L = 0.033;
D = 1.0;
m = 0.03;
x0 = zeros(3,1);
xN = [1 0 0].';

wall_pos = -0.01;

T = 5.0;
N = 20;
Ts = T/N;
    
load(['../../data_' num2str(Nm) '.mat'],'x0_init','xN_term');

% Number of variables
nx = (Nm-1)*2*3;
nu = 3;

EXPORT = 1;
COMPILE = 1;

ps = {}; vs = {}; states = [];
for i = 1:Nm-1
    DifferentialState(['p_' num2str(i) '(3)']);
    evalc(['p = p_' num2str(i)]);
    DifferentialState(['v_' num2str(i) '(3)']);
    evalc(['v = v_' num2str(i)]);
    ps = {ps{:}, p};
    vs = {vs{:}, v};
    states = [states; ps{i}];
    states = [states; vs{i}];
end
Control u(3);

% Compute forces
F = {};
for i = 1:Nm-1
    if i == 1
        dist = ps{1}-x0;
    else
        dist = ps{i}-ps{i-1};
    end
    tmp = D*(1 - L/sqrt(dist.'*dist));
    F = {F{:}, tmp*dist};
end

% Set up ODE
ode = [];
for i = 1:Nm-2
    f = 1/m*(F{i+1} - F{i}) - [0;0;g];
    ode = [ ode; dot(ps{i}) == vs{i} ];
    ode = [ ode; dot(vs{i}) == f ];
end
ode = [ ode; dot(ps{end}) == vs{end} ];
ode = [ ode; dot(vs{end}) == u ];

%% MPCexport
acadoSet('problemname', ['GN_' num2str(Nm)]);

ocp = acado.OCP( 0.0, T, N );

W_mat = eye(nx+nu,nx+nu);
WN_mat = eye(nx,nx);
W = acado.BMatrix(W_mat);
WN = acado.BMatrix(WN_mat);

ocp.minimizeLSQ( W, [states;u] );
ocp.minimizeLSQEndTerm( WN, states );

ocp.subjectTo( -10.0 <= u <= 10.0 );
for i = 1:Nm-1
    ocp.subjectTo( wall_pos <= ps{i}(2) );
end
ocp.subjectTo( 'AT_END', states == xN_term );

ocp.setModel(ode);

mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',          'FULL_CONDENSING_N2');
if d == 4
    method = '_GL8';
    mpc.set( 'INTEGRATOR_TYPE',         'INT_IRK_GL8'       );
elseif d == 3
    method = '_GL6';
    mpc.set( 'INTEGRATOR_TYPE',         'INT_IRK_GL6'       );
else
    error('These scripts were written for a 3- or 4-stage Gauss collocation method.')
end
mpc.set( 'IMPLICIT_INTEGRATOR_MODE', 	'IFT'               );
mpc.set( 'NUM_INTEGRATOR_STEPS',        Ns*N                );
mpc.set( 'QP_SOLVER',                   'QP_QPOASES'    	);
mpc.set( 'HOTSTART_QP',                 'NO'             	);

if EXPORT
    mpc.exportCode( ['export_GN_' num2str(Nm)] );
    copyfile('qpoases', ['export_GN_' num2str(Nm) '/qpoases'], 'f')
end
if COMPILE
    copyfile('rien_solver_mex.c', ['export_GN_' num2str(Nm) '/rien_solver_mex.c'], 'f')
    copyfile('make_rien_solver.m', ['export_GN_' num2str(Nm) '/make_rien_solver.m'], 'f')
    cd(['export_GN_' num2str(Nm)])
    make_rien_solver( ['../acado_GN_' num2str(Nm) '_' num2str(d)] )
    cd ..
end

%% PARAMETERS SIMULATION
input.x0 = x0_init.';
input.x = repmat(xN_term.',N+1,1);
input.u = zeros(N,nu);
input.od = [];
input.y = zeros(N,nx+nu);
input.yN = zeros(1,nx);

input.W = blkdiag(1e-10*eye(nx),eye(nu));
input.WN = blkdiag(1e-10*eye(nx));

%% SIMULATION LOOP
display('------------------------------------------------------------------')
display('               SQP Loop'                                    )
display('------------------------------------------------------------------')

cpuTime = []; simTime = []; qpTime = []; condenseTime = []; regularizeTime = [];

iter = 0;
KKT_val = 1;
outputs = {};
while KKT_val > 1e-12 && iter < 1000 
    
    % Solve NMPC OCP
    evalc(['output = acado_GN_' num2str(Nm) '_' num2str(d) '(input)']);
    outputs{iter+1} = output;
    
    cpuTime = [cpuTime output.info.cpuTime];
    simTime = [simTime output.info.simTime];
    qpTime = [qpTime output.info.qpTime];
    condenseTime = [condenseTime output.info.condenseTime];
    regularizeTime = [regularizeTime output.info.regularizeTime];
    
    % Save the MPC step
    input.x = output.x;
    input.u = output.u;
    
    KKT_val = output.info.kktValue;
    iter = iter+1;
    disp(['current iter: ' num2str(iter) '   ' char(9) ' (QP status: ' num2str(output.info.status) ', RTI step: ' num2str(output.info.cpuTime*1e6) ' µs, KKT val: ' num2str(KKT_val) ')'])
end

disp([' '])
disp(['Number of iterations   : ' num2str(iter)])
disp(['----------------------------------------------------'])
disp(['Average simTime        : ' num2str(round(mean(simTime)*1e6)*1e-3) 'ms'])
disp(['Average qpTime         : ' num2str(round(mean(qpTime)*1e6)*1e-3) 'ms'])
disp(['Average condenseTime   : ' num2str(round(mean(condenseTime)*1e6)*1e-3) 'ms'])
disp(['----------------------------------------------------'])
disp(['Average cpuTime        : ' num2str(round(mean(cpuTime)*1e6)*1e-3) 'ms'])

load(['../../data_ME_' num2str(Nm) method '.mat'],'x0_init','xN_term','resX','resU','res');

err_ACADO_x = max(max(abs(output.x-resX.')))
err_ACADO_u = max(max(abs(output.u-resU.')))

save(['../ACADO_GN_' num2str(Nm) method '.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime','outputs');

!rm *_RUN*
!rm *_data_acadodata*
end
