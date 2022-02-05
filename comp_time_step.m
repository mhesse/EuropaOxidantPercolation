function [dtDmax] = comp_time_step(phiD,qD,vD,XDc,ZDc,XDx,ZDx,XDz,ZDz,Grid,ve_cD)

% interpolate velocities to cell centers
% Fluid velocities
QDx = interp2(XDx,ZDx,reshape(qD(1:Grid.Nfx),Grid.Ny,Grid.Nx+1),XDc,ZDc);
QDz = interp2(XDz,ZDz,reshape(qD(Grid.Nfx+1:Grid.Nf),Grid.Ny+1,Grid.Nx),XDc,ZDc);

dtDxq = min(Grid.dx*phiD./(abs(QDx(:))));
dtDzq = min(Grid.dy*phiD./(abs(QDz(:))));
dtDmax.q = min([dtDxq,dtDzq]);
% Solid velocities
VDx = interp2(XDx,ZDx,reshape(vD(1:Grid.Nfx),Grid.Ny,Grid.Nx+1),XDc,ZDc);
VDz = interp2(XDz,ZDz,reshape(vD(Grid.Nfx+1:Grid.Nf),Grid.Ny+1,Grid.Nx),XDc,ZDc);
dtDxv = min(Grid.dx./(abs(VDx(:))));
dtDzv = min(Grid.dy./(abs(VDz(:))));
dtDmax.v = min([dtDxv,dtDzv]);

% Effective tracer velocity
Ve_cDx = interp2(XDx,ZDx,reshape(ve_cD(1:Grid.Nfx),Grid.Ny,Grid.Nx+1),XDc,ZDc);
Ve_cDz = interp2(XDz,ZDz,reshape(ve_cD(Grid.Nfx+1:Grid.Nf),Grid.Ny+1,Grid.Nx),XDc,ZDc);
dtDxve = min(Grid.dx./(abs(Ve_cDx(:))));
dtDzve = min(Grid.dy./(abs(Ve_cDz(:))));
dtDmax.ve = min([dtDxve,dtDzve]);

dtDmax.tot = min([dtDmax.v,dtDmax.q,dtDmax.ve]);

