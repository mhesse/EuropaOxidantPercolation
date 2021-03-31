function [hD,pD,qD] = solve_Helmholtz2D(D,G,I,phi,n,m,Grid,B,N,fn,Param,Yc)
% Porosity Matrices
Phi_n = comp_mean(reshape(phi.^n,Grid.Ny,Grid.Nx),-1,1,Grid);      % Note be careful expects a row vector!
Phi_m = spdiags(phi.^m,0,Grid.N,Grid.N);

%% Solve mod. Helmholtz equations 
L  = -D*Phi_n*G+Phi_m;              % system matrix
fs = Phi_m*Yc(:); % r.h.s.
flux = @(h) -Phi_n*G*h;
res = @(h,cell) L(cell,:)*h - fs(cell); 

%% Solve boundary value problem
hD = solve_lbvp(L,fs+fn,B,Param.g,N);
qD = comp_flux_gen(flux,res,hD,Grid,Param);
pD = hD-Yc(:);
end
