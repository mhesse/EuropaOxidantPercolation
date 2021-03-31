function [Ae] = comp_effective_flux_2nd(phiD,Aq,Av,dist_coeff,Grid,Scales)
% authors: Marc Hesse, Jake Jordan
% date: 13 Dec 2019
% Descritopn:
% Computes the effevtive velosity of a bulk quantity from solid velocity
% relative Darcy flux and the distribution coefficients.

% Input:
% phiD = scaled porosity
% Aq = advection matrix from the fluid flux
% Av = advection matrix from the solid velocity
% dist_coeff = distribution coefficient of the bulk quantity
phi_c = Scales.phi_c;
phiD_tilde = phiD+(1-phi_c*phiD).*dist_coeff/phi_c;
Ae = Aq*spdiags(1./phiD_tilde,0,Grid.N,Grid.N) + phi_c*Av;

    
