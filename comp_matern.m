function [m,gamma,delta] = comp_matern(sigma,rho,Grid,D,G,I,seed)

coeff = 3.6e-07; % handcrafted, locally sourced
% 6.3441e-06 - 1 by 1 domain
% 9.7629e-06 - 2 by 2 domain
% 9.0255e-06 - 3 by 3 domain
% 8.6570e-06 - 4 by 4 domain

gamma =sqrt(coeff)/2*rho/sigma;
delta = 2*sqrt(coeff)/(rho*sigma);

L = -gamma*D*G+delta*I;
rng(seed);
fs = randn(Grid.N,1);

Param.dof_dir = [];
Param.dof_f_dir = [];
Param.g = [];
Param.dof_neu = [];
Param.dof_f_neu = [];
Param.qb = [];

[B,N,fn] = build_bnd(Param,Grid,I);
m = solve_lbvp(L,fs+fn,B,Param.g,N);
end

