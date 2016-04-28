clc;
clear all;
close all;

import casadi.*

% Get collocation points
d = 4;
Ns = 3; % NUMBER OF INTEGRATION STEPS

for Nm = 3:8
    
    disp(['---- Nm value = ' num2str(Nm) '----']);
    
N = 20;
T_init = 10.0;
    
load(['../data_' num2str(Nm) '.mat'],'x0_init','xN_term');
    
% Environment
g = 9.81;     % [N/kg]
L = 0.033;
D = 1.0;
m = 0.03;
x_fixed = zeros(3,1);
xN = [1 0 0].';

wall_pos = -0.01;

T_min = 0.2;
T_max = 15;

% Number of variables
nx = (Nm-1)*2*3+1;
nu = 3;

% State variables
u = SX.sym('u',3);
dae.p = u;

dae.x = [];
states = [];
for i = 1:Nm-1
    p = SX.sym(['p' num2str(i)],3);
    v = SX.sym(['v' num2str(i)],3);
    

    x_struct = struct('p',p,'v',v);
    states = [states; x_struct];
    dae.x = [dae.x; casadi_struct2vec(x_struct)];
end
T = SX.sym('T');

% Compute forces
F = {};
for i = 1:Nm-1
    if i == 1
        dist = states(1).p-x_fixed;
    else
        dist = states(i).p-states(i-1).p;
    end
    tmp = D*(1 - L/sqrt(dist.'*dist));
    F = {F{:}, tmp*dist};
end

% Set up ODE
dae.ode = [];
for i = 1:Nm-2
    f = 1/m*(F{i+1} - F{i}) - [0;0;g];
    dae.ode = [dae.ode; casadi_vec(x_struct,'p',T*states(i).v,'v',T*f)];
end
dae.ode = [dae.ode; casadi_vec(x_struct,'p',T*states(end).v,'v',T*u)];

dae_full = dae;
dae_full.x = [dae_full.x; T];
dae_full.ode = [dae_full.ode; 0];

% Get collocation points
tau_root = collocationPoints(d,'legendre');

collfun = simpleColl(dae_full,tau_root,1/(Ns*N));
collfun = collfun.expand();

%% Find rest position
Xpoints = linspace(0,1,Nm);
x0_guess = [Xpoints(2:end);zeros(5,Nm-1)];
x0_guess = x0_guess(:);
u_guess = zeros(3,1);

odeFun = Function('odeFun',{dae.x,T,dae.p},{dae.ode,jacobian(dae.ode,[dae.x;dae.p])});
out = odeFun({x0_guess,1,u_guess});
while norm(full(out{1})) > 1e-10
    out = odeFun({x0_guess,1,u_guess});
    val = full(out{1});
    jac = full(out{2});
    delta = -jac\val;
    x0_guess = x0_guess + delta(1:nx-1);
end
% x0_guess
xN_term = x0_guess;
err_rest = norm(full(out{1}))

x0_mat2 = [zeros(6,1) reshape(x0_guess,6,Nm-1)];

Ypoints = linspace(0,1.5,Nm);
Zpoints = linspace(0,0.5,Nm);

x0_init = [zeros(1,Nm-1);Ypoints(2:end);Zpoints(2:end);zeros(3,Nm-1)];
x0_init = x0_init(:);
u_init = zeros(3,1);

out = odeFun({x0_init,1,u_init});
while norm(full(out{1})) > 1e-10
    out = odeFun({x0_init,1,u_init});
    val = full(out{1});
    jac = full(out{2});
    delta = -jac\val;
    x0_init = x0_init + delta(1:nx-1);
end
% x0_init
err_rest = norm(full(out{1}))

x0_mat = [zeros(6,1) reshape(x0_init,6,Nm-1)];

Fontsize = 20;
set(0,'DefaultAxesFontSize',Fontsize)

figure(1); set(gcf, 'Color','white');
plot3(x0_mat2(1,:), x0_mat2(2,:), x0_mat2(3,:), '--ro', 'MarkerSize', 10); hold on;
plot3(x0_mat(1,:), x0_mat(2,:), x0_mat(3,:), '--bo', 'MarkerSize', 10);
p = patch([0, 1, 1, 0], [wall_pos, wall_pos, wall_pos, wall_pos], [-4, -4, 1, 1], 'g');
xlabel( 'x [m]');
ylabel( 'y [m]');
zlabel( 'z [m]');
xlim([0 1]);
ylim([-0.1 2]);
zlim([-4 1]);
title('Initial and reference point');
legend('reference', 'initial','wall')
view([145 25]);
grid on;

%% Optimal Control Problem
Xs = {};
for i=1:N+1
   Xs{i} = MX.sym(['X_' num2str(i)],nx);
end
XCs = {};
Us = {};
for i=1:N
   XCs{i} = MX.sym(['XC_' num2str(i)],nx,Ns*d);
   Us{i}  = MX.sym(['U_' num2str(i)],nu);
end

V_block = struct();
V_block.X  = Sparsity.dense(nx,1);
V_block.XC  = Sparsity.dense(nx,Ns*d);
V_block.U  = Sparsity.dense(nu,1);

% Simple bounds on states
lbx = {};
ubx = {};

% List of constraints
g = {};

% List of all decision variables (determines ordering)
V = {};
for k=1:N
  % Add decision variables
  V = {V{:} casadi_vec(V_block,'X',Xs{k},'XC',XCs{k},'U',Us{k})};
  
  if k==1
    % Bounds at t=0
    x_lb = [x0_init;T_min];
    x_ub = [x0_init;T_max];
    u_lb = -10*ones(3,1);
    u_ub = 10*ones(3,1);
    lbx = {lbx{:} casadi_vec(V_block,-inf,'X',x_lb,'U',u_lb)};
    ubx = {ubx{:} casadi_vec(V_block,inf, 'X',x_ub,'U',u_ub)};
  else %k < N
    % Bounds for other t
    m_lb  = [-inf;wall_pos;-inf;-inf;-inf;-inf];
    m_ub  = [inf;inf;inf;inf;inf;inf];
    x_lb = [repmat(m_lb,Nm-1,1);T_min];
    x_ub = [repmat(m_ub,Nm-1,1);T_max];
    u_lb = -10.0*ones(3,1);
    u_ub = 10.0*ones(3,1);
    lbx = {lbx{:} casadi_vec(V_block,-inf,'X',x_lb,'U',u_lb)};
    ubx = {ubx{:} casadi_vec(V_block,inf, 'X',x_ub,'U',u_ub)};
  end
  
  % Obtain collocation expressions
  Xcur = Xs{k};
  for i = 1:Ns
    coll_out = collfun({Xcur,XCs{k}(:,1+(i-1)*d:i*d),Us{k}});
    Xcur = coll_out{1};
    g = {g{:} coll_out{2}};         % collocation constraints
  end
  g = {g{:} Xs{k+1}-Xcur}; % gap closing
end
  
V = {V{:} Xs{end}};

% Bounds for final t
x_lb = [xN_term;T_min];
x_ub = [xN_term;T_max];
lbx = {lbx{:} x_lb};
ubx = {ubx{:} x_ub};

% Objective function
controls = vertcat(Us{:});

objective = Xs{1}(end);

nlp = struct('x',vertcat(V{:}), 'f', objective, 'g', vertcat(g{:}));

nlpfun = Function('nlp',nlp,char('x','p'),char('f','g'));

opts.ipopt = struct('linear_solver','ma27','acceptable_tol', 1e-12, 'tol', 1e-12);
% opts.ipopt = struct('acceptable_tol', 1e-10, 'tol', 1e-10);
solver = nlpsol('solver','ipopt',nlp, opts);

u_guess = zeros(3,1);
x0 = [repmat([xN_term;T_init;repmat([xN_term;T_init],Ns*d,1);u_guess],N,1);xN_term;T_init];

args = struct;
args.x0 = x0;
args.lbx = vertcat(lbx{:});
args.ubx = vertcat(ubx{:});
args.lbg = 0;
args.ubg = 0;

res = solver(args);

x0 = full(res.x);

dim = size(casadi_struct2vec(V_block));
res_split = vertsplit(res.x,dim(1));

res_U = {};
for r=res_split(1:end-1)
    rs = casadi_vec2struct(V_block,r{1});
    res_U = {res_U{:} rs.U};
end

res_X = {};
for r=res_split(1:end-1)
    rs = casadi_vec2struct(V_block,r{1});
    res_X = {res_X{:} rs.X};
end
res_X = {res_X{:} res_split{end}};

% Visualization solution
T_opt = full(res_X{1}(end))

figure(2); set(gcf, 'Color','white');
plot3(x0_mat2(1,:), x0_mat2(2,:), x0_mat2(3,:), '--ro', 'MarkerSize', 14); hold on;
plot3(x0_mat(1,:), x0_mat(2,:), x0_mat(3,:), '--bo', 'MarkerSize', 14);
p = patch([0, 1, 1, 0], [wall_pos, wall_pos, wall_pos, wall_pos], [-4, -4, 1, 1], 'g');
for k = 1:N+1
    tmpX = full(res_X{k});
    x_mat = [zeros(6,1) reshape(tmpX(1:end-1),6,Nm-1)];
    plot3(x_mat(1,:), x_mat(2,:), x_mat(3,:), ':k+', 'MarkerSize', 6);
end
xlabel( 'x [m]');
ylabel( 'y [m]');
zlabel( 'z [m]');
xlim([0 1]);
ylim([-0.1 2]);
zlim([-4 1]);
title('Initial and reference point');
legend('reference','initial','wall','solution')
view([145 25]);
grid on;

resX = vertcat(res_X{:});
resX = full(reshape(resX,nx,N+1));

resU = vertcat(res_U{:});
resU = full(reshape(resU,nu,N));

res = full(res.x);

end
