clc;
clear all;
close all;

Ns = 3;
d = 4;
method = '';

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

T_init = 2.0;
T_min = 0.2;
T_max = 15;

wall_pos = -0.01;

T = 1.0;
N = 20;
Ts = T/N;
    
load(['../../data_' num2str(Nm) '.mat'],'x0_init','xN_term');

% Number of variables
nx = (Nm-1)*2*3+1;
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
DifferentialState T_total;
Control u(3);

% Compute forces
F = {};
for i = 1:Nm-1
    if i == 1
        dist = is(ps{1}-x0);
    else
        dist = is(ps{i}-ps{i-1});
    end
    tmp = is(D*(1 - L/sqrt(dist.'*dist)));
    F = {F{:}, is(tmp*dist)};
end

% Set up ODE
ode = [];
for i = 1:Nm-2
    f = 1/m*(F{i+1} - F{i}) - [0;0;g];
    ode = [ ode; dot(ps{i}) == T_total*vs{i} ];
    ode = [ ode; dot(vs{i}) == T_total*f ];
end
ode = [ ode; dot(ps{end}) == T_total*vs{end} ];
ode = [ ode; dot(vs{end}) == T_total*u ];
ode = [ ode; dot(T_total) == 0 ];

%% MPCexport
acadoSet('problemname', ['EH_EN_' num2str(Nm)]);

ocp = acado.OCP( 0.0, T, N );

ocp.minimizeLagrangeTerm( T_total );

ocp.subjectTo( -10.0 <= u <= 10.0 );
for i = 1:Nm-1
    ocp.subjectTo( wall_pos <= ps{i}(2) );
end
ocp.subjectTo( T_min <= T_total <= T_max );
ocp.subjectTo( 'AT_START', states == x0_init );
ocp.subjectTo( 'AT_END', states == xN_term );

ocp.setModel(ode);

mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',       'EXACT_HESSIAN'     );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',          'CONDENSING_N2'     );  % <--- TIME OPTIMAL FORMULATION
if d == 4
    method = '_GL8';
    mpc.set( 'INTEGRATOR_TYPE',         'INT_IRK_GL8'       );
elseif d == 3
    method = '_GL6';
    mpc.set( 'INTEGRATOR_TYPE',         'INT_IRK_GL6'       );
else
    error('These scripts were written for a 3- or 4-stage Gauss collocation method.')
end
mpc.set( 'IMPLICIT_INTEGRATOR_MODE', 	'LIFTED' 			);
mpc.set( 'LIFTED_INTEGRATOR_MODE',      1                   );
mpc.set( 'DYNAMIC_SENSITIVITY',         'SYMMETRIC'         );
mpc.set( 'NUM_INTEGRATOR_STEPS',        Ns*N                );
mpc.set( 'QP_SOLVER',                   'QP_QPOASES'    	);
mpc.set( 'HOTSTART_QP',                 'NO'             	);
mpc.set( 'FIX_INITIAL_STATE',           'NO'				);   % <--- TIME OPTIMAL FORMULATION

if EXPORT
    mpc.exportCode( ['export_EH_EN_' num2str(Nm)] );
    copyfile('qpoases', ['export_EH_EN_' num2str(Nm) '/qpoases'], 'f')
end
if COMPILE
    copyfile('rien_solver_mex_EN.c', ['export_EH_EN_' num2str(Nm) '/rien_solver_mex.c'], 'f')
    copyfile('make_rien_solver.m', ['export_EH_EN_' num2str(Nm) '/make_rien_solver.m'], 'f')
    cd(['export_EH_EN_' num2str(Nm)])
    make_rien_solver( ['../acado_EH_EN_' num2str(Nm) '_' num2str(d)] )
    cd ..
end

%% PARAMETERS SIMULATION
load(['../ACADO_GN_' num2str(Nm) method '.mat'],'input');

input.x0 = [];
input.W = [];
input.WN = [];
input.x = [input.x repmat(T_init,N+1,1)];
input.y = zeros(N,nx+nu);
input.yN = zeros(1,nx);

input.mu = zeros(N,nx);

input.control = 1; % INITIALIZE
evalc(['output = acado_EH_EN_' num2str(Nm) '_' num2str(d) '(input)']);
input.control = 0; % NORMAL BEHAVIOUR

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
    evalc(['output = acado_EH_EN_' num2str(Nm) '_' num2str(d) '(input)']);
    outputs{iter+1} = output;
    
    cpuTime = [cpuTime output.info.cpuTime];
    simTime = [simTime output.info.simTime];
    qpTime = [qpTime output.info.qpTime];
    condenseTime = [condenseTime output.info.condenseTime];
    regularizeTime = [regularizeTime output.info.regularizeTime];
    
    % Save the MPC step
    input.x = output.x;
    input.u = output.u;
    input.mu = output.mu;
    
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
disp(['Average regularizeTime : ' num2str(round(mean(regularizeTime)*1e6)*1e-3) 'ms'])
disp(['----------------------------------------------------'])
disp(['Average cpuTime        : ' num2str(round(mean(cpuTime)*1e6)*1e-3) 'ms'])

load(['../../data_TO_' num2str(Nm) method '.mat'],'x0_init','xN_term','resX','resU','res');

err_ACADO_x = max(max(abs(output.x-resX.')))
err_ACADO_u = max(max(abs(output.u-resU.')))

save(['../ACADO_EH_EN_' num2str(Nm) method '.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime','outputs');

!rm *_RUN*
!rm *_data_acadodata*
end
