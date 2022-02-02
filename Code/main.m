%%% Approximately compute the 2nd-order energy correction term (containing all blocks)
%%% using mode-wise JL embeddings. Comparison of Gaussian, RFD, RCD, and
%%% Rademacher JL. A 2nd-stage JL can be performed if chosen to do so.

clear

addpath('./Raw_Data')

% =============================================================
% inputs
% =============================================================
Interaction='em1.8-2.0'; % 'em1.8-2.0', 'em2.0-2.0', 'em2.0-2.5', 'em2.8-2.0', ...
atom='Sn132';
name_suff='_hwHO024';
eMax=6;
lMax=10; % Matters for eMax>10

dim_frac=[0.1 0.15 0.2 0.25 0.3 0.4 0.5];

N_trials=50; % number of trials

JL_mats={'Gaussian','RFD'}; % 1st-stage JL matrices. A set from {'Gaussian','RFD','RCD','Rademacher'}
secJL_flag=1; % set to 0 if secondary JL should NOT be done
dim_frac_sec=0.2;

sec_JL_type='RFD'; % 'Gaussian' or 'RFD'

% =============================================================
% =============================================================

name_pref=['hint_' Interaction '_'];
if eMax>lMax
    eMax_lMax_str=sprintf('_eMax%02d_lMax%02d', eMax, lMax);
else
    eMax_lMax_str=sprintf('_eMax%02d', eMax);
end

J_max=2*min(lMax, eMax)+1;
n_J=J_max+1;
suff='hf';
disp('Loading denominator...')
load(['D_' name_pref atom suff eMax_lMax_str name_suff '.mat'])

t_start=tic;

siz=size(D);
n_dims=length(siz);

if any(strcmp(JL_mats, 'Gaussian'))
    E2_proj=zeros(N_trials, length(dim_frac), n_J);
end
if any(strcmp(JL_mats, 'RFD'))
    E2_proj_RFD=zeros(N_trials, length(dim_frac), n_J);
end
if any(strcmp(JL_mats, 'RCD'))
    E2_proj_RCD=zeros(N_trials, length(dim_frac), n_J);
end
if any(strcmp(JL_mats, 'Rademacher'))
    E2_proj_Rad=zeros(N_trials, length(dim_frac), n_J);
end
if secJL_flag~=0
    if any(strcmp(JL_mats, 'Gaussian'))
        E2_proj_sec=zeros(N_trials, length(dim_frac), n_J);
    end
    if any(strcmp(JL_mats, 'RFD'))
        E2_proj_RFD_sec=zeros(N_trials, length(dim_frac), n_J);
    end
    if any(strcmp(JL_mats, 'RCD'))
        E2_proj_RCD_sec=zeros(N_trials, length(dim_frac), n_J);
    end
    if any(strcmp(JL_mats, 'Rademacher'))
        E2_proj_Rad_sec=zeros(N_trials, length(dim_frac), n_J);
    end
    d_sec=zeros(1, length(dim_frac));
end

E2_proj_tot=cell(1, length(JL_mats));
for J=0:n_J-1
    
    fprintf('Channel number %d \t t_elapsed=%f\n', J, toc(t_start))
    
    H=load([name_pref atom suff eMax_lMax_str '_J'  num2str(J,'%02d') name_suff '.mat']);
    H=H.data;
    
    H_perm=permute(H,[3 4 1 2]);
    
        
    %%% Calculate modewise random projections
    d=zeros(1, n_dims);
    for p=1:length(dim_frac)
        for kk=1:n_dims
            d(kk)=ceil(dim_frac(p)*siz(kk));
        end
        if secJL_flag~=0
            d_sec(p)=ceil(dim_frac_sec*prod(d));
            R_rand_idx_sec=floor(prod(d)*rand(1,d_sec(p)));
        end
        for k=1:N_trials            
            fprintf('\nJ=%d \t dim_frac(%d) \t Trial %d \n', J, p, k)
            toc(t_start)
            
            %%% construct projection matrices
            if any(strcmp(JL_mats, 'Gaussian'))
                A=cell(1,n_dims);
            end
            if any(strcmp(JL_mats, 'RFD'))
                A_RFD=cell(1,n_dims);
            end
            if any(strcmp(JL_mats, 'RCD'))
                A_RCD=cell(1,n_dims);
            end
            if any(strcmp(JL_mats, 'Rademacher'))
                A_Rad=cell(1,n_dims);
            end
            for kk=1:n_dims
                if any(strcmp(JL_mats, 'Gaussian'))
                    A{kk}=randn(d(kk), siz(kk));
                    A{kk}=A{kk}/sqrt(d(kk)); % normalize the columns of A
                end
                
                R=datasample(1:siz(kk), d(kk)); % random restriction with replacement
                D_rademacher=rademacher(1, siz(kk));
                if any(strcmp(JL_mats, 'RFD'))
                    A_RFD{kk}=FastJLmat_RFD(R, siz(kk), D_rademacher);
                end
                if any(strcmp(JL_mats, 'RCD'))
                    A_RCD{kk}=FastJLmat_RCD(R, siz(kk), D_rademacher);
                end
                if any(strcmp(JL_mats, 'Rademacher'))
                    A_Rad{kk}=rademacher(d(kk), siz(kk))/sqrt(d(kk));
                end
            end
            
            if secJL_flag~=0
                if strcmp(sec_JL_type, 'Gaussian')
                    A_sec=randn(d_sec(p), prod(d))/sqrt(d_sec(p));
                elseif strcmp(sec_JL_type, 'RFD')
                    D_sec=rademacher(1, prod(d));
                end
            end
            
            %%% compute the approximate Energy terms
            
            % Gaussian JL
            if any(strcmp(JL_mats, 'Gaussian'))
                H1=tensor_randproj_multi_opt(A, 'full', H);
                H2=tensor_randproj_multi_opt(A, 'full', H_perm.*D);
                E2_proj(k,p,J+1)=-0.25*(2*J+1)*H1(:)'*H2(:);
                E2_proj_all{1}(k,p,J+1)=E2_proj(k,p,J+1);
            
                if secJL_flag~=0
                    if strcmp(sec_JL_type, 'Gaussian')
                        E2_proj_sec(k,p,J+1)=-0.25*(2*J+1)*(A_sec*H1(:))'*(A_sec*H2(:));
                    elseif strcmp(sec_JL_type, 'RFD')
                        H1=H1(:);
                        H2=H2(:);
                        H1=fft(sparse(1:length(H1), 1:length(H1), D_sec)*H1);
                        H1=H1(R_rand_idx_sec+1)/sqrt(d_sec(p));
                        H2=fft(sparse(1:length(H2), 1:length(H2), D_sec)*H2);
                        H2=H2(R_rand_idx_sec+1)/sqrt(d_sec(p));
                        
                        E2_proj_sec(k,p,J+1)=-0.25*(2*J+1)*H1'*H2;
                        E2_proj_sec_all{1}(k,p,J+1)=E2_proj_sec(k,p,J+1);
                    end
                end
            end
            
            % RFD JL
            if any(strcmp(JL_mats, 'RFD'))
                H1=tensor_randproj_multi_opt(A_RFD, 'full', H);
                H2=tensor_randproj_multi_opt(A_RFD, 'full', H_perm.*D);
                E2_proj_RFD(k,p,J+1)=-0.25*(2*J+1)*H1(:)'*H2(:);
                E2_proj_all{2}(k,p,J+1)=E2_proj_RFD(k,p,J+1);
            
                if secJL_flag~=0
                    if strcmp(sec_JL_type, 'Gaussian')
                        E2_proj_RFD_sec(k,p,J+1)=-0.25*(2*J+1)*(A_sec*H1(:))'*(A_sec*H2(:));
                    elseif strcmp(sec_JL_type, 'RFD')
                        H1=H1(:);
                        H2=H2(:);
                        H1=fft(sparse(1:length(H1), 1:length(H1), D_sec)*H1);
                        H1=H1(R_rand_idx_sec+1)/sqrt(d_sec(p));
                        H2=fft(sparse(1:length(H2), 1:length(H2), D_sec)*H2);
                        H2=H2(R_rand_idx_sec+1)/sqrt(d_sec(p));
                        
                        E2_proj_RFD_sec(k,p,J+1)=-0.25*(2*J+1)*H1'*H2;
                        E2_proj_sec_all{2}(k,p,J+1)=E2_proj_RFD_sec(k,p,J+1);
                    end
                end
            end
            
            % RCD JL
            if any(strcmp(JL_mats, 'RCD'))
                H1=tensor_randproj_multi_opt(A_RCD, 'full', H);
                H2=tensor_randproj_multi_opt(A_RCD, 'full', H_perm.*D);
                E2_proj_RCD(k,p,J+1)=-0.25*(2*J+1)*H1(:)'*H2(:);
                E2_proj_all{3}(k,p,J+1)=E2_proj_RCD(k,p,J+1);
                
                if secJL_flag~=0
                    if strcmp(sec_JL_type, 'Gaussian')
                        E2_proj_RCD_sec(k,p,J+1)=-0.25*(2*J+1)*(A_sec*H1(:))'*(A_sec*H2(:));
                    elseif strcmp(sec_JL_type, 'RFD')
                        H1=H1(:);
                        H2=H2(:);
                        H1=fft(sparse(1:length(H1), 1:length(H1), D_sec)*H1);
                        H1=H1(R_rand_idx_sec+1)/sqrt(d_sec(p));
                        H2=fft(sparse(1:length(H2), 1:length(H2), D_sec)*H2);
                        H2=H2(R_rand_idx_sec+1)/sqrt(d_sec(p));
                        
                        E2_proj_RCD_sec(k,p,J+1)=-0.25*(2*J+1)*H1'*H2;
                        E2_proj_sec_all{3}(k,p,J+1)=E2_proj_RCD_sec(k,p,J+1);
                    end
                end
            end
            
            % Rademacher
            if any(strcmp(JL_mats, 'Rademacher'))
                H1=tensor_randproj_multi_opt(A_Rad, 'full', H);
                H2=tensor_randproj_multi_opt(A_Rad, 'full', H_perm.*D);
                E2_proj_Rad(k,p,J+1)=-0.25*(2*J+1)*H1(:)'*H2(:);
                E2_proj_all{4}(k,p,J+1)=E2_proj_Rad(k,p,J+1);
                
                if secJL_flag~=0
                    if strcmp(sec_JL_type, 'Gaussian')
                        E2_proj_Rad_sec(k,p,J+1)=-0.25*(2*J+1)*(A_sec*H1(:))'*(A_sec*H2(:));
                    elseif strcmp(sec_JL_type, 'RFD')
                        H1=H1(:);
                        H2=H2(:);
                        H1=fft(sparse(1:length(H1), 1:length(H1), D_sec)*H1);
                        H1=H1(R_rand_idx_sec+1)/sqrt(d_sec(p));
                        H2=fft(sparse(1:length(H2), 1:length(H2), D_sec)*H2);
                        H2=H2(R_rand_idx_sec+1)/sqrt(d_sec(p));
                        
                        E2_proj_Rad_sec(k,p,J+1)=-0.25*(2*J+1)*H1'*H2;
                        E2_proj_sec_all{4}(k,p,J+1)=E2_proj_Rad_sec(k,p,J+1);
                    end
                end
            end
                           
            
        end
    end
    
end

end_time=toc(t_start);

E2_proj_all(cellfun(@isempty,E2_proj_all))=[];
if secJL_flag~=0
    E2_proj_sec_all(cellfun(@isempty,E2_proj_sec_all))=[];
end

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

%%% save results

if secJL_flag~=0
    save(['results_JL_E2_sum_sec_' Interaction '_' atom suff eMax_lMax_str name_suff,...
        sprintf('_%s_compare4_Ntrials%d_dimfrac(1)%.2d_dimfrac(end)%.2d_dimfracsec%.2d.mat',...
        sec_JL_type, N_trials, dim_frac(1)*1e3, dim_frac(end)*1e3, dim_frac_sec*1e3)],...
        'E2_proj_all', 'E2_proj_sec_all',...
        'dim_frac', 'dim_frac_sec', 'N_trials',...
        'end_time', 'eMax', 'secJL_flag', 'sec_JL_type',...
        'JL_mats', 'atom', 'Interaction')
else
    save(['results_JL_E2_sum_' Interaction '_' atom suff eMax_lMax_str name_suff,...
        sprintf('_compare4_Ntrials%d_dimfrac(1)%.2d_dimfrac(end)%.2d.mat',...
        N_trials, dim_frac(1)*1e3, dim_frac(end)*1e3)],...
        'E2_proj_all',...
        'dim_frac', 'N_trials', 'end_time', 'eMax',...
        'secJL_flag', 'JL_mats', 'atom', 'Interaction')
end

