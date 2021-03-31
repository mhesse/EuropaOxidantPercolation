function [CD] = evolve_concentration_total(D,I,Av,CD,phiD,phiD_old,qD,B,N,Param,Grid,Scales,theta,dtD,Dist)
% author: Marc Hesse
% date: 17 Apr 2019
phi_c = Scales.phi_c;

% Update concentration
phi_tilde = @(phi) phi + (1-phi_c*phi)*Dist/phi_c;
% phi_tilde = @(phi) phi_c*phi + (1-phi_c*phi)*Dist;
phi_theta = theta*phi_tilde(phiD_old) + (1-theta)*phi_tilde(phiD);
Aq = flux_upwind(qD,Grid);
L_C = D*(Aq*spdiags(1./phi_theta,0,Grid.N,Grid.N)+phi_c*Av);
IM_C = I + dtD*(1-theta)*L_C;
EX_C = I - dtD*theta    *L_C;
CD = solve_lbvp(IM_C,EX_C*CD,B,Param.g,N);