function [theta] = smoothness(df,q,Grid)
% author: Marc Hesse
% date: 29 Jul 2016
% Description: Computes an indicator of the smoothness of the solution at
%              cell interfaces from divided differences. For details see
%              LeVeque (2002) Chp. 6.11 on p. 114.
%
% Input: 
% df = Nf by 1 vector of difference of cell centered values on faces
% q  = Nf by 1 vector of fluxes on faces
% Grid = structure containing all pertinent information about the grid.
%
% Output: 
% theta = Nf by 1 vector of divided difference in upwind direction
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I] = build_ops(Grid);
% >> c = sin(Grid.xc);
% >> q = ones(Grid.Nfx,1);
% >> dc = Grid.dx*G*c;
% >> theta = smoothness(dc,q,Grid);

Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;
Nfx = Grid.Nfx; % # of x faces
Nfy = Grid.Nfy; % # of y faces
Nf = length(q);
dx = Grid.dx;

tol = 1e-10;
high_value = 100;
%% Compute divided differences
if ((Nx>1) && (Ny==1)) || ((Nx==1) && (Ny>1)) % 1D
    Nf = length(q); %N = length(df)
    theta = zeros(Nf,1);
    
    % Upwing divided difference
    df_up = (q>=0).*[0;df(1:Nf-1)] + (q<0).*[df(2:Nf);0];
    
    %% Deal with problematic cases: df << 1 denominator blows up
    tol = 1e-10;
    
    % Case 1: Both df, and df_up small => function is smooth => theta = 1
    ids = find(abs(df)+abs(df_up)<tol); % both jumps small
    theta(ids) = 1;
    
    % Case 2: % is df small and df_up relatively large
    ids = find(abs(df)<tol & abs(df_up)>tol );
    theta(ids) = high_value; % choose theta large
    
    % Case 3: df is finite
    ids = find(abs(df)>= tol);
    theta(ids) = df_up(ids)./df(ids);
else % 2D
    %% Smoothness in x-direction
    theta_x = zeros(Nfx,1);
    DX = reshape(df(1:Nfx),    Ny  ,Nx+1);
    DX_up = (reshape(q(1:Grid.Nfx),Ny,Nx+1)>=0).*[zeros(Ny,1),DX(:,1:Nx)] +...
            (reshape(q(1:Grid.Nfx),Ny,Nx+1)< 0).*[DX(:,2:Nx+1),zeros(Ny,1)];
    
    % Case 1: Both dfx, and dfx_up small => function is smooth => theta = 1
    ids = find(abs(df(1:Nfx))+abs(DX_up(:))<tol);
    theta_x(ids) = 1; clear ids
    
    % Case 2: % is dfx small and dfx_up relatively large => choose theta large
    ids = find(abs(df(1:Nfx))<tol & abs(DX_up(:))>tol );
    theta_x(ids) = high_value; clear ids
    
    % Case 3: df is finite
    ids = find(abs(df(1:Nfx))>= tol);
    theta_x(ids) = DX_up(ids)./df(ids); clear ids % Note: This works because x-differences are first in the df vector!
    
    %% Smoothness in the y-direction
    theta_y = zeros(Nfy,1);
    DY = reshape(df(Nfx+1:Nf),Ny+1,Nx  );
    DY_up = (reshape(q(Nfx+1:Nf),Ny+1,Nx)>=0).*[zeros(1,Nx);DY(1:Ny,:)] +...
            (reshape(q(Nfx+1:Nf),Ny+1,Nx)< 0).*[DY(2:Ny+1,:);zeros(1,Nx)];
        
    % Case 1: Both dfy, and dfy_up small => function is smooth => theta = 1
    ids = find(abs(df(Nfx+1:Nf))+abs(DY_up(:))<tol);
    theta_y(ids) = 1; clear ids
    
    % Case 2: % is dfy small and dfy_up relatively large => choose theta large
    ids = find(abs(df(Nfx+1:Nf))<tol & abs(DY_up(:))>tol );
    theta_y(ids) = high_value; clear ids
    
    % Case 3: df is finite
    ids = find(abs(df(Nfx+1:Nf))>= tol);
    theta_y(ids) = DY_up(ids)./df(Nfx+ids); clear ids % Note: Need to shift index because the x-differences are first in the df vector!

    theta = [theta_x;theta_y];
end

if sum(isnan(theta))+ sum(isinf(theta)) > 0 
    error('Smoothness indicator for flux limiter is not computed correctly.')

end
