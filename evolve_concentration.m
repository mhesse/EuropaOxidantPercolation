function [cD,CD] = evolve_concentration(D,Av,cD,phiD,phiD_old,qD,B,N,Param,Grid,Scales,theta,dtD,Dist)
% author: Marc Hesse
% date: 17 Apr 2019
phi_c = Scales.phi_c;

% Update concentration
phi_tilde = @(phi) phi + (1-phi_c*phi)*Dist/phi_c;
phi_theta = theta*phi_tilde(phiD_old) + (1-theta)*phi_tilde(phiD);
Aq = flux_upwind(qD,Grid);
L = D*(Aq+phi_c*Av*spdiags(phi_theta,0,Grid.N,Grid.N));
IM = spdiags(phiD,    0,Grid.N,Grid.N) + dtD*(1-theta)*L;
EX = spdiags(phiD_old,0,Grid.N,Grid.N) - dtD*theta    *L;
cD = solve_lbvp(IM,EX*cD,B,Param.g,N);
CD = phiD.*cD;