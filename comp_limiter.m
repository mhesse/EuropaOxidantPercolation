function [phi_lim] = comp_limiter(theta,limiter)
% author: Marc Hesse
% date: 29 Jul 2016
% Description: Computes the slope limiter for 2nd-order TVD scheme.  
%              For details see LeVeque (2002) Chp. 6.11 on p. 115.
% 
% Input: 
% theta = Nf by 1 vector of divided differences 
% limiter = indicates type of slope limiter:
%           gd = Godunov 
%           lw = Lax-Wendroff
%           bw = Beam Warming
%           fr = Fromm
%           mm = minmod
%           sb = suberbee
%           mc = MC-limiter (monotenized central-difference limiter)
%           vl = van Leer
%           va = van Alba

% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I] = build_ops(Grid);
% >> c = sin(Grid.xc);
% >> q = ones(Grid.Nfx,1);
% >> dc = Grid.dx*G*c;
% >> theta = smoothness(dc,q,Grid);
% >> phi_lim = comp_limiter(theta,'mc');

if strcmp(limiter,'gd') % Godunov/Upwind
    phi_lim = spalloc(length(theta),1,0) ; 
elseif strcmp(limiter,'lw') % Lax-Wendroff
    phi_lim = ones(length(theta),1);
elseif strcmp(limiter,'bw') % Beam-Warming
    phi_lim = theta;
elseif strcmp(limiter,'fr') % Fromm
    phi_lim = .5*(1+theta);
elseif strcmp(limiter,'mm') % minmod
    phi_lim = minmod(1,theta);
elseif strcmp(limiter,'sb') % superbee
    phi_lim = max(0,max(min(1,2*theta),min(2,theta)));
elseif strcmp(limiter,'mc') % monotonized central difference
    phi_lim = max(0,min(min((1+theta)/2,2*theta),2));
elseif strcmp(limiter,'vl') % van Leer
    phi_lim = ( theta + abs(theta) )./( 1 + abs(theta) );
elseif strcmp(limiter,'va') % van Alba
    phi_lim = ( theta.^2 + theta )./( theta.^2 + 1 );
    phi_lim(theta < 0) = 0;
else
    error('Unknown slope limiter.')
end
    