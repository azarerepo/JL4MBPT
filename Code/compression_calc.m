function [comp_1st, comp_2nd]=compression_calc(dim_frac, dim_frac_sec, eMax, implicit_numdim)
%%% implicit_dim: implicit number of dimensions of tensors forming the
%%% inner product

eMax_vals=[2,4,6,8,10,12,14];
mode_sizes=[12,30,56,90,132,174,216];
dim_size=mode_sizes(eMax_vals==eMax);

target_dim=ceil(dim_frac*dim_size).^implicit_numdim; % no. of dimensions  of implicit tensors in the inner product
comp_1st=target_dim/dim_size^implicit_numdim; % implicit tensors in the inner product are 5D
target_dim_sec=ceil(target_dim*dim_frac_sec);
comp_2nd=target_dim_sec/dim_size^implicit_numdim;