clc;
clear all;
close all;

T = 5; N = 20; Ts = T/N;

Fontsize = 22;
set(0,'DefaultAxesFontSize',Fontsize)
MarkerSize = 14;

Nm = 5;

%% minimum control effort
load(['data_ME_' num2str(Nm) '_GL6.mat'],'x0_init','xN_term','resX','resU','res');
res_x_ME_GL6 = resX; res_u_ME_GL6 = resU;
res_x_ME = resX; res_u_ME = resU;


load(['ACADO_SQP/ACADO_GN_' num2str(Nm) '_GL6.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime','outputs');
errs_GN = [];
for i = 1:length(outputs)
    err_x = max(max(abs(outputs{i}.x-res_x_ME_GL6.'))); err_u = max(max(abs(outputs{i}.u-res_u_ME_GL6.')));
    errs_GN = [errs_GN; max(err_x,err_u)];
end

load(['ACADO_SQP/ACADO_GN_EN_' num2str(Nm) '_GL6.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime','outputs');
errs_GN_EN = [];
for i = 1:length(outputs)
    err_x = max(max(abs(outputs{i}.x-res_x_ME_GL6.'))); err_u = max(max(abs(outputs{i}.u-res_u_ME_GL6.')));
    errs_GN_EN = [errs_GN_EN; max(err_x,err_u)];
end

load(['ACADO_SQP/ACADO_GN_IN_' num2str(Nm) '_GL6.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime','outputs');
errs_GN_IN = [];
for i = 1:length(outputs)
    err_x = max(max(abs(outputs{i}.x-res_x_ME_GL6.'))); err_u = max(max(abs(outputs{i}.u-res_u_ME_GL6.')));
    errs_GN_IN = [errs_GN_IN; max(err_x,err_u)];
end

load(['ACADO_SQP/ACADO_GN_INIS_' num2str(Nm) '_GL6.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime','outputs');
errs_GN_INIS = [];
for i = 1:length(outputs)
    err_x = max(max(abs(outputs{i}.x-res_x_ME_GL6.'))); err_u = max(max(abs(outputs{i}.u-res_u_ME_GL6.')));
    errs_GN_INIS = [errs_GN_INIS; max(err_x,err_u)];
end

TOL = 2e-7;

%% time optimal
load(['data_TO_' num2str(Nm) '_GL6.mat'],'x0_init','xN_term','resX','resU','res');
res_x_TO = resX(1:end-1,:); res_u_TO = resU;
T2 = resX(end,1); Ts2 = T2/N;


load(['ACADO_SQP/ACADO_EH_' num2str(Nm) '_GL6.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime','outputs');
errs_EH = [];
for i = 1:length(outputs)
    err_x = max(max(abs(outputs{i}.x-resX.'))); err_u = max(max(abs(outputs{i}.u-resU.')));
    errs_EH = [errs_EH; max(err_x,err_u)];
end


load(['ACADO_SQP/ACADO_EH_EN_' num2str(Nm) '_GL6.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime','outputs');
errs_EH_EN = [];
for i = 1:length(outputs)
    err_x = max(max(abs(outputs{i}.x-resX.'))); err_u = max(max(abs(outputs{i}.u-resU.')));
    errs_EH_EN = [errs_EH_EN; max(err_x,err_u)];
end


load(['ACADO_SQP/ACADO_EH_INIS_' num2str(Nm) '_GL6.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime','outputs');
errs_EH_INIS = [];
for i = 1:length(outputs)
    err_x = max(max(abs(outputs{i}.x-resX.'))); err_u = max(max(abs(outputs{i}.u-resU.')));
    errs_EH_INIS = [errs_EH_INIS; max(err_x,err_u)];
end


load(['ACADO_SQP/ACADO_EH_IN_' num2str(Nm) '_GL6.mat'],'input','output','iter','simTime','qpTime','condenseTime','regularizeTime','cpuTime','outputs');
errs_EH_IN = [];
for i = 1:length(outputs)
    err_x = max(max(abs(outputs{i}.x-resX.'))); err_u = max(max(abs(outputs{i}.u-resU.')));
    errs_EH_IN = [errs_EH_IN; max(err_x,err_u)];
end


%% figures
figure;
pos = {'p_x End', 'p_y End', 'p_z End'};
ref = xN_term(end-5:end-3);
for k = 1:3
    subplot(3,2,2*(k-1)+1);
    plot([0:Ts:T], res_x_ME(end-6+k,:), '--bx', 'MarkerSize', 10); hold on;
    plot([0:Ts2:T2], res_x_TO(end-6+k,:), '--ro', 'MarkerSize', 10);
    plot([0 T], ref(k)*ones(1,2), '--k');
    xlabel('time (s)')
    ylabel(pos(k))
    if k == 3
    legend('control effort', 'time optimal')
    end
end
us = {'u_x', 'u_y', 'u_z'};
for k = 1:3
    subplot(3,2,2*k);
    stairs([0:Ts:T-Ts], res_u_ME(k,:), '--bx', 'MarkerSize', 10); hold on;
    stairs([0:Ts2:T2-Ts2], res_u_TO(k,:), '--ro', 'MarkerSize', 10);
    xlabel('time (s)')
    ylabel(us(k))
    if k == 3
    legend('control effort', 'time optimal')
    end
end

figure;
subplot(2,1,1);
semilogy(errs_GN, '--kx', 'MarkerSize', MarkerSize); hold on;
semilogy(errs_GN_EN, '--ro', 'MarkerSize', MarkerSize);
semilogy(errs_GN_IN, '--gs', 'MarkerSize', MarkerSize);
semilogy(errs_GN_INIS(1:end-3), '--b+', 'MarkerSize', MarkerSize);
legend('Without lifting (MS)', 'Exact lifting (LC-EN)', 'IN lifting (LC-IN)', 'INIS lifting (LC-AF-INIS)')
title('Minimum effort: Gauss-Newton based SQP')
ylim([1e-10 1e1])
xlim([1 50])
xlabel('iteration number')
ylabel('||W-W^*||_\infty')

subplot(2,1,2);
semilogy(errs_EH(1:end-2), '--kx', 'MarkerSize', MarkerSize); hold on;
semilogy(errs_EH_EN(1:end-2), '--ro', 'MarkerSize', MarkerSize);
semilogy(errs_EH_IN(1:end-5), '--gs', 'MarkerSize', MarkerSize);
semilogy(errs_EH_INIS(1:end-3), '--b+', 'MarkerSize', MarkerSize);
legend('Without lifting (MS)', 'Exact lifting (LC-EN)', 'IN lifting (LC-IN)', 'INIS lifting (LC-INIS)')
title('Time optimal: Exact Hessian based SQP')
ylim([1e-8 1e2])
xlim([1 20])
xlabel('iteration number')
ylabel('||W-W^*||_\infty')


