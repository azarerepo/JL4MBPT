
clear

% load the saved results and true energy correction values.
load('results_JL_E2_sum_sec_em1.8-2.0_Sn132hf_eMax06_hwHO024_RFD_compare4_Ntrials50_dimfrac(1)100_dimfrac(end)500_dimfracsec200_1.mat')
load('E2_true_em1.8-2.0_Sn132_eMax06_hwHO024.mat')

lMax=10; % plays a role when eMax > 10

ise = evalin( 'base', 'exist(''E2_proj_all'',''var'') == 1' );
if ise
    for k=1:length(JL_mats)
        JL_mat=JL_mats{k};
        if strcmp(JL_mat, 'Gaussian')
            E2_proj=E2_proj_all{k};
            if secJL_flag~=0
                E2_proj_sec=E2_proj_sec_all{k};
            end
        elseif strcmp(JL_mat, 'RFD')
            E2_proj_RFD=E2_proj_all{k};
            if secJL_flag~=0
                E2_proj_RFD_sec=E2_proj_sec_all{k};
            end
        elseif strcmp(JL_mat, 'RCD')
            E2_proj_RCD=E2_proj_all{k};
            if secJL_flag~=0
                E2_proj_RCD_sec=E2_proj_sec_all{k};
            end
        elseif strcmp(JL_mat, 'Rad')
            E2_proj_Rad=E2_proj_all{k};
            if secJL_flag~=0
                E2_proj_Rad_sec=E2_proj_sec_all{k};
            end
        end
    end
    clear E2_proj_all E2_proj_sec_all
end

E2_orig1=sum(E2_orig);

J_max=2*min(lMax, eMax)+1;
n_J=J_max+1;

[comp_1st, comp_2nd]=compression_calc(dim_frac, dim_frac_sec, eMax, 4);

%%%%% Plots below include Gaussian and RFD embeddings. The same type of
%%%%% plots can be made when RCD and Rademacher JL embeddings are used.

%%%%==================================================================================
%%%% Gaussian JL
%%%%==================================================================================

%%% plots of relative error with mean and standard deviation and box plots
E2_proj2=sum(E2_proj,3);
Err_rel2=(E2_proj2-E2_orig1*ones(N_trials, length(dim_frac)))./(E2_orig1*ones(N_trials, length(dim_frac)));

%%%%==================================================================================
%%%% RFD JL
%%%%==================================================================================

E2_proj1_RFD=sum(E2_proj_RFD,3);
E2_proj1_RFD=mean(E2_proj1_RFD);
E2_rel_RFD=E2_proj1_RFD./E2_orig1;


%%% plots of relative error with mean and standard deviation and box plots
E2_proj2_RFD=sum(E2_proj_RFD,3);
Err_rel2_RFD=(E2_proj2_RFD-E2_orig1*ones(N_trials, length(dim_frac)))./(E2_orig1*ones(N_trials, length(dim_frac)));


%%%%=================================================================================
%%%% plots for comparing the Gaussian and RFD
%%%%=================================================================================

%%% Mean absolute
figure,
plot(comp_1st, mean(abs((sum(E2_proj,3)-E2_orig1)./E2_orig1)), '*', 'MarkerSize', 12, 'LineWidth', 1.4)
hold on
plot(comp_1st, mean(abs(real((sum(E2_proj_RFD,3)-E2_orig1)./E2_orig1))), '*', 'MarkerSize', 12, 'LineWidth', 1.4) % absolute value of the real part
xlabel('$c_{tot}$','Interpreter','latex', 'FontSize', 14)
ylabel('$\overline{\Delta E^{(2)}}$','Interpreter','latex', 'FontSize', 14)
% ylh = get(gca,'ylabel');
% set(ylh, 'Rotation', 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
title(sprintf('%s, $eMax=%d$ (%d trials).',atom, eMax, N_trials),'Interpreter','latex', 'FontSize', 14)
grid on
legend('Gaussian','RFD_{real}', 'FontSize', 14)
set(gca, 'XScale', 'Log')


%%% only 2nd JL
figure,
plot(comp_2nd, mean(abs((sum(E2_proj_sec,3)-E2_orig1)./E2_orig1)), '*', 'MarkerSize', 12, 'LineWidth', 1.4)
hold on
plot(comp_2nd, mean(abs(real((sum(E2_proj_RFD_sec,3)-E2_orig1)./E2_orig1))), '*', 'MarkerSize', 12, 'LineWidth', 1.4) % absolute value of the real part
xlabel('$c_{tot}$','Interpreter','latex', 'FontSize', 14)
ylabel('$\overline{\Delta E^{(2)}}$','Interpreter','latex', 'FontSize', 14)
% ylh = get(gca,'ylabel');
% set(ylh, 'Rotation', 0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
title(sprintf('%s, 2nd JL, $eMax=%d$ (%d trials).', atom, eMax, N_trials),'Interpreter','latex', 'FontSize', 14)
grid on
legend(['(Gaussian + '  sprintf('%s)_{real}', sec_JL_type)],['(RFD + ' sprintf('%s)_{real}', sec_JL_type)], 'FontSize', 14)
set(gca, 'XScale', 'Log')

%%% 2nd JL and modewise JL
figure,
plot(comp_2nd, mean(abs((sum(E2_proj_sec,3)-E2_orig1)./E2_orig1)), 'd', 'MarkerSize', 12, 'LineWidth', 1.4)
hold on
plot(comp_2nd, mean(abs(real((sum(E2_proj_RFD_sec,3)-E2_orig1)./E2_orig1))), 'd', 'MarkerSize', 12, 'LineWidth', 1.4) % absolute value of the real part
hold on
plot(comp_1st, mean(abs((sum(E2_proj,3)-E2_orig1)./E2_orig1)), '*', 'MarkerSize', 12, 'LineWidth', 1.4)
hold on
plot(comp_1st, mean(abs(real((sum(E2_proj_RFD,3)-E2_orig1)./E2_orig1))), '*', 'MarkerSize', 12, 'LineWidth', 1.4) % absolute value of the real part
xlabel('$c_{tot}$','Interpreter','latex', 'FontSize', 14)
ylabel('$\overline{\Delta E^{(2)}}$','Interpreter','latex', 'FontSize', 14)
title(sprintf('%s, $c_2=%.2f, eMax=%d$ (%d trials).', atom, dim_frac_sec, eMax, N_trials),'Interpreter','latex', 'FontSize', 14)
grid on
legend(['(Gaussian + '  sprintf('%s)_{real}', sec_JL_type)],['(RFD + '  sprintf('%s)_{real}', sec_JL_type)],'Gaussian','RFD_{real}', 'FontSize', 14)
set(gca, 'XScale', 'Log')

figure,
plot(comp_2nd, mean(abs(real((sum(E2_proj_RFD_sec,3)-E2_orig1)./E2_orig1))), 'd', 'MarkerSize', 12, 'LineWidth', 1.4) % absolute value of the real part
hold on
plot(comp_1st, mean(abs(real((sum(E2_proj_RFD,3)-E2_orig1)./E2_orig1))), '*', 'MarkerSize', 12, 'LineWidth', 1.4) % absolute value of the real part
xlabel('$c_{tot}$','Interpreter','latex', 'FontSize', 14)
ylabel('$\overline{\Delta E^{(2)}}$','Interpreter','latex', 'FontSize', 14)
title(sprintf('%s, $c_2=%.2f, eMax=%d$ (%d trials).', atom, dim_frac_sec, eMax, N_trials),'Interpreter','latex', 'FontSize', 14)
grid on
legend(['(RFD + '  sprintf('%s)_{real}', sec_JL_type)],'RFD_{real}', 'FontSize', 14)
set(gca, 'XScale', 'Log')
%%%%=================================================================================
%%%% box plots for selected compression ratios
%%%%=================================================================================
idx_sel=[length(dim_frac)-3 length(dim_frac)];%[4 7];
comp_1st_sel=comp_1st(idx_sel);
comp_2nd_sel=comp_2nd(idx_sel);
y_min=-2; y_max=2;


%%% box plots of absolute errors (real part for RFD) with both compressions on the horizontal axis
figure,
boxplot([abs(Err_rel2(:,idx_sel(1))), abs(real(Err_rel2_RFD(:,idx_sel(1)))),...
    abs(Err_rel2(:,idx_sel(2))), abs(real(Err_rel2_RFD(:,idx_sel(2))))],...
    'Notch','on','Labels',...
    {sprintf('$c_{tot}=%.5f$',comp_1st_sel(1)),...
    sprintf('$c_{tot}=%.5f$',comp_1st_sel(1)),...
    sprintf('$c_{tot}=%.5f$',comp_1st_sel(2)),...
    sprintf('$c_{tot}=%.5f$',comp_1st_sel(2))}, 'colors', 'brbr')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 13)
ylabel('$\left| \left(E_p^{(2)}-E^{(2)}\right)/E^{(2)} \right|$','Interpreter','latex', 'FontSize', 14)
title(sprintf('%s, eMax=%d (%d trials).', atom, eMax, N_trials),'Interpreter','latex', 'FontSize', 14)
boxes=findobj(gca,'Tag','Box');
legend(boxes([2 1]),'Gaussian','RFD_{real}', 'FontSize', 14)

%%% secondary JL: box plots of absolute errors (real part for RFD)
dim_frac_sel=dim_frac(idx_sel);
E2_proj2_sec=sum(E2_proj_sec,3);
Err_rel2_sec=(E2_proj2_sec-E2_orig1*ones(N_trials, length(dim_frac)))./(E2_orig1*ones(N_trials, length(dim_frac)));
E2_proj2_RFD_sec=sum(E2_proj_RFD_sec,3);
Err_rel2_RFD_sec=(E2_proj2_RFD_sec-E2_orig1*ones(N_trials, length(dim_frac)))./(E2_orig1*ones(N_trials, length(dim_frac)));
figure,
boxplot([abs(Err_rel2_sec(:,idx_sel(1))), abs(real(Err_rel2_RFD_sec(:,idx_sel(1)))),...
    abs(Err_rel2_sec(:,idx_sel(2))), abs(real(Err_rel2_RFD_sec(:,idx_sel(2))))],...
    'Notch','on','Labels',...
    {sprintf('$c_{tot}=%.5f$',comp_2nd_sel(1)),...
    sprintf('$c_{tot}=%.5f$',comp_2nd_sel(1)),...
    sprintf('$c_{tot}=%.5f$',comp_2nd_sel(2)),...
    sprintf('$c_{tot}=%.5f$',comp_2nd_sel(2))}, 'colors', 'brbr')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 13)
ylabel('$\left| \left(E_p^{(2)}-E^{(2)}\right)/E^{(2)} \right|$','Interpreter','latex', 'FontSize', 14)
title(sprintf('%s, $2^{nd}$ JL, $c_2$=%.2f, eMax=%d (%d trials).', atom, dim_frac_sec, eMax, N_trials),'Interpreter','latex', 'FontSize', 14)
boxes=findobj(gca,'Tag','Box');
legend(boxes([2 1]),['(Gaussian + '  sprintf('%s)_{real}', sec_JL_type)],['(RFD + '  sprintf('%s)_{real}', sec_JL_type)], 'FontSize', 14)


figure,
boxplot([abs(real(Err_rel2_RFD_sec(:,idx_sel(1)))),abs(real(Err_rel2_RFD_sec(:,idx_sel(2))))],...
    'Notch','on','Labels',...
    {sprintf('$c_{tot}=%.5f$',comp_2nd_sel(1)),...
    sprintf('$c_{tot}=%.5f$',comp_2nd_sel(2))})
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 13)
ylabel('$\left| \left(E_p^{(2)}-E^{(2)}\right)/E^{(2)} \right|$','Interpreter','latex', 'FontSize', 14)
title(sprintf('%s, $2^{nd}$ JL, $c_2$=%.2f, eMax=%d (%d trials).', atom, dim_frac_sec, eMax, N_trials),'Interpreter','latex', 'FontSize', 14)



