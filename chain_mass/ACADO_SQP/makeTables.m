clc;
clear all;
close all;

method = '_GL8';
nm = 5;

%% GAUSS-NEWTON
disp('------    GAUSS-NEWTON AVERAGE COMPUTATION TIMES (in ms)    ------')
disp('                   standard    lifting    IN-adjoint    INIS-single  ')
err_mat = []; sim_mat = []; qp_mat = []; condense_mat = []; cpu_mat = [];
for Nm = 3:7
    
    errs = []; simTimes = []; qpTimes = []; condenseTimes = []; cpuTimes = [];
    
    load(['../data_ME_' num2str(Nm) method '.mat'],'x0_init','xN_term','resX','resU','res');
    
    load(['ACADO_GN_' num2str(Nm) method '.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime');
    err_x = max(max(abs(output.x-resX.')));
    err_u = max(max(abs(output.u-resU.')));
    errs = [errs max(err_x,err_u)];
    simTimes = [simTimes mean(simTime)];
    qpTimes = [qpTimes mean(qpTime)];
    condenseTimes = [condenseTimes mean(condenseTime)];
    cpuTimes = [cpuTimes mean(cpuTime)];
    
    load(['ACADO_GN_EN_' num2str(Nm) method '.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime');
    err_x = max(max(abs(output.x-resX.')));
    err_u = max(max(abs(output.u-resU.')));
    errs = [errs max(err_x,err_u)];
    simTimes = [simTimes mean(simTime)];
    qpTimes = [qpTimes mean(qpTime)];
    condenseTimes = [condenseTimes mean(condenseTime)];
    cpuTimes = [cpuTimes mean(cpuTime)];
    
    load(['ACADO_GN_IN_' num2str(Nm) method '.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime');
    err_x = max(max(abs(output.x-resX.')));
    err_u = max(max(abs(output.u-resU.')));
    errs = [errs max(err_x,err_u)];
    simTimes = [simTimes mean(simTime)];
    qpTimes = [qpTimes mean(qpTime)];
    condenseTimes = [condenseTimes mean(condenseTime)];
    cpuTimes = [cpuTimes mean(cpuTime)];
    
    load(['ACADO_GN_INIS_' num2str(Nm) method '.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime');
    err_x = max(max(abs(output.x-resX.')));
    err_u = max(max(abs(output.u-resU.')));
    errs = [errs max(err_x,err_u)];
    simTimes = [simTimes mean(simTime)];
    qpTimes = [qpTimes mean(qpTime)];
    condenseTimes = [condenseTimes mean(condenseTime)];
    cpuTimes = [cpuTimes mean(cpuTime)];
    
    err_mat = [err_mat;errs];
    sim_mat = [sim_mat;simTimes];
    qp_mat = [qp_mat;qpTimes];
    condense_mat = [condense_mat;condenseTimes];
    cpu_mat = [cpu_mat;cpuTimes];
    
    nx = (Nm-1)*2*3;
    disp([' Nm = ' num2str(Nm) ', nx = ' num2str(nx) ':   ' num2str(round(cpuTimes*1e5)*1e-2)])
end
% err_mat

disp(' ')
disp(' ')
disp(['------    GAUSS-NEWTON DETAILED TIMINGS Nm = ' num2str(nm) ' (in ms)    ------'])
disp('                 standard   lifting   IN-adjoint   INIS-single  ')
disp(['simulation      ' num2str(round(sim_mat(nm-2,:)*1e5)*1e-2)])
disp(['qp condensing   ' num2str(round(condense_mat(nm-2,:)*1e5)*1e-2)])
disp(['qp solution     ' num2str(round(qp_mat(nm-2,:)*1e5)*1e-2)])
disp('---------------------------------------------------------')
disp(['total SQP step  ' num2str(round(cpu_mat(nm-2,:)*1e5)*1e-2)])

disp(' ')
disp(' ')


%% EXACT-HESSIAN
disp('------    EXACT-HESSIAN AVERAGE COMPUTATION TIMES (in ms)    ------')
disp('                   standard    lifting    IN-adjoint    INIS-single  ')
err_mat = []; sim_mat = []; qp_mat = []; condense_mat = []; reg_mat = []; cpu_mat = [];
for Nm = 3:7
    
    errs = []; simTimes = []; qpTimes = []; condenseTimes = []; regularizeTimes = []; cpuTimes = [];
    
    load(['../data_TO_' num2str(Nm) method '.mat'],'x0_init','xN_term','resX','resU','res');
    
    load(['ACADO_EH_' num2str(Nm) method '.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime');
    err_x = max(max(abs(output.x-resX.')));
    err_u = max(max(abs(output.u-resU.')));
    errs = [errs max(err_x,err_u)];
    simTimes = [simTimes mean(simTime)];
    qpTimes = [qpTimes mean(qpTime)];
    condenseTimes = [condenseTimes mean(condenseTime)];
    regularizeTimes = [regularizeTimes mean(regularizeTime)];
    cpuTimes = [cpuTimes mean(cpuTime)];
    
    load(['ACADO_EH_EN_' num2str(Nm) method '.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime');
    err_x = max(max(abs(output.x-resX.')));
    err_u = max(max(abs(output.u-resU.')));
    errs = [errs max(err_x,err_u)];
    simTimes = [simTimes mean(simTime)];
    qpTimes = [qpTimes mean(qpTime)];
    condenseTimes = [condenseTimes mean(condenseTime)];
    regularizeTimes = [regularizeTimes mean(regularizeTime)];
    cpuTimes = [cpuTimes mean(cpuTime)];
    
    load(['ACADO_EH_IN_' num2str(Nm) method '.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime');
    err_x = max(max(abs(output.x-resX.')));
    err_u = max(max(abs(output.u-resU.')));
    errs = [errs max(err_x,err_u)];
    simTimes = [simTimes mean(simTime)];
    qpTimes = [qpTimes mean(qpTime)];
    condenseTimes = [condenseTimes mean(condenseTime)];
    regularizeTimes = [regularizeTimes mean(regularizeTime)];
    cpuTimes = [cpuTimes mean(cpuTime)];
    
    load(['ACADO_EH_INIS_' num2str(Nm) method '.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime');
    err_x = max(max(abs(output.x-resX.')));
    err_u = max(max(abs(output.u-resU.')));
    errs = [errs max(err_x,err_u)];
    simTimes = [simTimes mean(simTime)];
    qpTimes = [qpTimes mean(qpTime)];
    condenseTimes = [condenseTimes mean(condenseTime)];
    regularizeTimes = [regularizeTimes mean(regularizeTime)];
    cpuTimes = [cpuTimes mean(cpuTime)];
    
    err_mat = [err_mat;errs];
    sim_mat = [sim_mat;simTimes];
    qp_mat = [qp_mat;qpTimes];
    condense_mat = [condense_mat;condenseTimes];
    reg_mat = [reg_mat;regularizeTimes];
    cpu_mat = [cpu_mat;cpuTimes];
    
    nx = (Nm-1)*2*3+1;
    disp([' Nm = ' num2str(Nm) ', nx = ' num2str(nx) ':   ' num2str(round(cpuTimes*1e5)*1e-2)])
end
% err_mat

disp(' ')
disp(' ')
disp(['------    Exact-Hessian DETAILED TIMINGS Nm = ' num2str(nm) ' (in ms)    ------'])
disp('                 standard   lifting   IN-adjoint   INIS-single  ')
disp(['simulation      ' num2str(round(sim_mat(nm-2,:)*1e5)*1e-2)])
disp(['qp condensing   ' num2str(round(condense_mat(nm-2,:)*1e5)*1e-2)])
disp(['regularization  ' num2str(round(reg_mat(nm-2,:)*1e5)*1e-2)])
disp(['qp solution     ' num2str(round(qp_mat(nm-2,:)*1e5)*1e-2)])
disp('---------------------------------------------------------')
disp(['total SQP step  ' num2str(round(cpu_mat(nm-2,:)*1e5)*1e-2)])

disp(' ')
disp(' ')
