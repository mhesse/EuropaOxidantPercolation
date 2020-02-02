% CH_OxidantTransport2D_limiters.m
% author: Marc Hesse, Jake Jordan
% date: 30 Jan 2020
% Description:
% This code solves the Europa oxidant transport problem (high pororsity and tracer at top) with heat
% conduction and phase behavior in 2D. This is the last step before this
% was CH_OxidantTransport2D.m.
close all, clear, clc

%% Define physical parameters
% set_demo_defaults

% Domain Geometry
Param.refine = 2;
Param.dim.H = 25e3;            % Ice shell thickness  [m]
Param.dim.L = 10e3;            % Width of the domain  [m]
Param.dim.D = 0e3;             % Depth of ocean below [m]
Param.dim.melt.depth = 1.5e3;  % Center of melt lens  [m]
Param.dim.melt.thick  = 3e3;   % Thickness of melt lens  [m]
Param.exp.m = 1;               % compaction viscosity exponent [-]
Param.conc1.Dist = 0e-1;       % Dist = 0 -> all in fluid
% Param.perm = 'vonBargen';
% Param.perm = 'WarkWatson';
Param.perm = 'GHP';

Param.dim.z_tbl = 1e3;         % Base of thermal boundary layer [m]
Param.dim.z_sal = 5e3;         % Base of salt enriched crust [m]

% Initial condition
ic.z_tbl  = 24.9e3;       % Base of the thermal boundary layer  [m]
ic.z_sal  = 22e3;         % Base of the enriched crust [m]
ic.X_top = 0.01/0.263;    % Salt (HH) mass fraction in enriched crust [1]
ic.X_con = 1e-4/0.263;    % Salt mass fraction in convecting shell [1]
ic.Ttop = 100;            % Mean surface temperature [K]
ic.Tbot = 262;            % Basal temperature [K]
ic.seed = 2;              % seed for random perturbations

% Thermodynamic properties
Phase.L  = 3.34e5;
Phase.Xe = 0.8859; % for Hydrohalite
Phase.T1 = 273.15;
Phase.Te = 252.05;
[rho,cp,kappa,hD,HD_from_TPhi,Phase] = physical_properties_eutectic(Phase);
[HD_from_TX,CD_from_TX,phi_of_TX,f_of_TX,Phase] = setup_phase_behavior(rho,cp,hD,Phase);

% Characteristic scales of the ice shell
TD_c = (ic.Tbot-Phase.Te)/Phase.DT;
Scales.d_c   = 1e-3;  % grain size [m]
Scales.T_c   = 262;   % average temperature in conv. shell [K]

[TD,Phi,reg] = eval_phase_behavior(HD_from_TX(TD_c,ic.X_con),CD_from_TX(TD_c,ic.X_con),rho,cp,hD,Phase);
Scales.phi_c = Phi(3);  % background porosity [-]

[Scales,Param,ConstFuns] = comp_dimless_params_oxidant(Scales,Param);
[Scales,Param,TD_fun] = comp_dimless_params_enthalpy(Scales,Param,Phase,rho,cp,kappa);

%% Dimensionless parameters
Param.non_dim.tmax = 1.5;
TD_top = TD_fun(ic.Ttop);
TD_bot = TD_fun(ic.Tbot);

%% Define Grid and build operators
Grid.xmin = 0;  Grid.xmax = Param.non_dim.L; Grid.Nx = Param.refine*ceil(Param.non_dim.L);
Grid.ymin = -Param.non_dim.D;  Grid.ymax = Param.non_dim.H; Grid.Ny = Param.refine*ceil(Param.non_dim.H+Param.non_dim.D);
Nt = Param.non_dim.tmax*6e2;
Param.plot_interval = 5;

theta = 1; % forward  Euler (explicit) - Note: going implicit screws up conservation

Grid    = build_grid(Grid);
[D,G,I] = build_ops(Grid);
[XDc,ZDc] = meshgrid(Grid.xc,Grid.yc);
dtD     = Param.non_dim.tmax/Nt;
% zDc = Grid.yc; zDf = Grid.yf;

%% Define initial condition
[CD,HD,TD,Phi,phiD,CD_tracer,cD_tracer,X] = make_ic_oxidant_2D(ic,Scales,Grid,Param,ConstFuns,TD_fun,ZDc,XDc,HD_from_TX,CD_from_TX,rho,cp,hD,Phase,D,G,I);
dist_tracer = Param.conc1.Dist*ones(Grid.N,1);  

%% Define boundary conditions
% 1) Helmholtz eqn.
Param_h.dof_dir   = [];
Param_h.dof_f_dir = [];
Param_h.g         = [];
Param_h.dof_neu   = [];
Param_h.dof_f_neu = [];
Param_h.qb        = [];
[B_h,N_h,fn_h] = build_bnd(Param_h,Grid,I);

% 2) Poisson eqn.
Param_u.dof_dir   = [1];
Param_u.dof_f_dir = [1];
Param_u.g         = [0];
Param_u.dof_neu   = [];
Param_u.dof_f_neu = [];
Param_u.qb        = [];
[B_u,N_u,fn_u] = build_bnd(Param_u,Grid,I);

% 3) Advection eqn. - Bulk Concentration
Param_C.dof_dir   = [];
Param_C.dof_f_dir = [];
Param_C.g         = [];
Param_C.dof_neu   = [];
Param_C.dof_f_neu = [];
Param_C.qb        = [];
[B_C,N_C,fn_C] = build_bnd(Param_C,Grid,I);

% 4) Advection eqn. - Bulk Enthalpy
Param_H.dof_dir   = [Grid.dof_ymin;Grid.dof_ymax];
Param_H.dof_f_dir = [Grid.dof_f_ymin;Grid.dof_f_ymax];
Param_H.g = [HD_from_TPhi(TD_bot,Phi(Grid.dof_ymin,:));...
             HD_from_TPhi(TD_top,Phi(Grid.dof_ymax,:))];
Param_H.dof_neu   = [];
Param_H.dof_f_neu = [];
Param_H.qb        = [];
[B_H,N_H,fn_H] = build_bnd(Param_H,Grid,I);

% 5) Advection eqn. - Oxidant
Param_c.dof_dir   = [];
Param_c.dof_f_dir = [];
Param_c.g         = [];
Param_c.dof_neu   = [];
Param_c.dof_f_neu = [];
Param_c.qb        = [];
Param_c.dof_out = [];
[B_c,N_c,fn_c] = build_bnd(Param_c,Grid,I);

%% Initial condition
tD = 0; i = 0; k = 0;
phiD(abs(phiD)<1e-12) = 1e-5; 
[headD,pD,qD] = solve_Helmholtz(D,G,I,phiD,Param.exp.n,Param.exp.m,Grid,B_h,N_h,fn_h,Param_h,ZDc);
[uD,vD]    = solve_Poisson(D,G,I,phiD,Param.exp.m,pD,Grid,B_u,N_u,fn_u,Param_u);
figure('units','normalized','outerposition',[0 0 .5 1])
plot_solution_CH_2D(Grid,phiD,cD_tracer,CD,HD,TD,CD_tracer,hD,uD,qD,vD,pD,ZDc,XDc,Scales,Phase)

save(['solution_0.mat'],'tD','i','k','phiD','CD','headD','uD','qD','vD','pD','CD_tracer','cD_tracer','TD','Phi','Grid','Scales','Phase','ic','XDc','ZDc')

% Time evolution variables
phi_total = zeros(Nt+1,1); phi_total(1) = sum(phiD*Grid.dy);
Ox_total   = zeros(Nt+1,1);   Ox_total(1) = sum(CD_tracer*Grid.dy);
Ox_ocean = zeros(Nt+1,1);
time = zeros(Nt+1,1);
Pe = Param.non_dim.Pe;
dtD_ref = dtD; 

while tD < Param.non_dim.tmax
    %% find time step
    tD = tD + dtD; %dtD*i;
%     dt_max_v = max([Grid.dx./vD(1:Grid.Nfx);Grid.dy./vD(Grid.Nfx+1:Grid.Nf)])
%     dt_max_v = max([Grid.dx./vD(1:Grid.Nfx);Grid.dy./vD(Grid.Nfx+1:Grid.Nf)])
    
    %% Phase behavior and distribution coefficients
    [TD,Phi,reg] = eval_phase_behavior(HD,CD,rho,cp,hD,Phase);
    [dist_CD,dist_HD,dist_tracer] = comp_dist_coeffs(CD,HD,TD,Phi,dist_tracer,reg,rho,hD,Phase);
    % Temperature/Enthalpy BC's
    TD(Param_H.dof_dir) = [TD_bot*ones(Grid.Nx,1);TD_top*ones(Grid.Nx,1)];
    HD(Param_H.dof_dir) = [HD_from_TPhi(TD_bot*ones(Grid.Nx,1),Phi(Grid.dof_ymin,:));...
                           HD_from_TPhi(TD_top*ones(Grid.Nx,1),Phi(Grid.dof_ymax,:))];
    
    % Move tracer out of frozen region
%     top = find(Phi(:,3),1,'last'); % top of melting region
%     CD_tracer(top) = CD_tracer(top)+sum(CD_tracer(Phi(:,3)==0));
%     CD_tracer(Phi(:,3)==0) = 0;                   
                       
    % Compute total mass of fluid and tracer
    phiD = Phi(:,3)/Scales.phi_c;
    phi_total(i+1) = sum(phiD*Grid.dx*Grid.dy); C_total(i+1) =  sum(CD*Grid.dx*Grid.dy); time(i+1) = tD;
    
    % Compute tracer in ocean
    Ox_ocean(i+1) =  sum(CD_tracer(ZDc(:)<0)*Grid.dx*Grid.dy)/Ox_total(1)*100;

    %% Flow calculations
    phiD(abs(phiD)<1e-12) = 1e-5; 
    [headD,pD,qD] = solve_Helmholtz(D,G,I,phiD,Param.exp.n,Param.exp.m,Grid,B_h,N_h,fn_h,Param_h,ZDc);
    [uD,vD]    = solve_Poisson(D,G,I,phiD,Param.exp.m,pD,Grid,B_u,N_u,fn_u,Param_u);
    
    %% Transport calculations
    % Update porosity
    Aq = flux_upwind(qD,Grid);
    Av = flux_upwind(vD,Grid);
    Aq_tracer = build_adv_op(qD,CD_tracer,dtD,G,Grid,Param_c,'mc');
    Av_tracer = build_adv_op(vD,CD_tracer,dtD,G,Grid,Param_c,'mc');
    
    Ae_CD = comp_effective_flux(phiD,Aq,Av,dist_CD,Grid,Scales);
    Ae_HD = comp_effective_flux(phiD,Aq,Av,dist_HD,Grid,Scales);
    Ae_cD = comp_effective_flux(phiD,Aq_tracer,Av_tracer,dist_tracer,Grid,Scales);
    
    
    Ex_CD = I - dtD*D*Ae_CD;
    Ex_cD = I - dtD*D*(Ae_cD);
    CD = solve_lbvp(I,Ex_CD*CD,B_C,Param_C.g,N_C);
    CD_tracer = solve_lbvp(I,Ex_cD*CD_tracer,B_c,Param_c.g,N_c);
    cD_tracer = CD_tracer./ConstFuns.phi_tilde(phiD);
    
    % Direct update of HD in mixed HD-TD formulation
    el_kap = reshape(kappa.sysD(TD,Phi),Grid.Ny,Grid.Nx);
    Kappa = comp_mean(el_kap,-1,1,Grid);
    HD = HD - dtD*D*(Ae_HD*HD - 1/Param.non_dim.PeT*Kappa*G*TD);    
    
    if mod(i,Param.plot_interval) == 0
        k = k+1;
        fprintf('total: i = %d: tD= %3.2f;\n',i,tD);
        plot_solution_CH_2D(Grid,phiD,cD_tracer,CD,HD,TD,CD_tracer,hD,uD,qD,vD,pD,ZDc,XDc,Scales,Phase)
        save(['solution_',num2str(k),'.mat'],'tD','i','k','phiD','CD','hD','uD','qD','vD','pD','CD_tracer','cD_tracer','TD','Phi')
    end
end

%% Save time evolution
filename = ['CHTransport_','_Ny',num2str(Grid.Ny),'Nt',num2str(Nt),'_m',num2str(Param.exp.m),'.mat'];
Ny = Grid.Ny; m = Param.exp.m;
save(filename,'time','phi_total','C_total','Nt','Ny','m','Param','theta','phiD','Grid','Ox_ocean')

%% Plot oxidant arrival in ocean
figure
plot(time*Scales.t_c/(60^2*24*365.25),Ox_ocean)
xlabel 'time in years', ylabel 'Oxidant in ocean [%]'


function plot_solution_CH_2D(Grid,phiD,cD,CD,HD,TD,CD_tracer,hD,uD,qD,vD,pD,ZDc,XDc,Scales,Phase)
T = Phase.Te + (Phase.T1-Phase.Te)*TD;

subplot 121
surf(XDc,ZDc,Scales.phi_c*reshape(phiD,Grid.Ny,Grid.Nx)), view(2), shading interp, colorbar
xlabel 'x_D'
ylabel 'z_D'
title '\phi_D'
axis equal tight

subplot 122
surf(XDc,ZDc,reshape(CD_tracer,Grid.Ny,Grid.Nx)), view(2), shading interp, colorbar
xlabel 'x_D'
ylabel 'z_D'
title 'C_D tracer'
axis equal tight

% plot_HC_diagram(Phase,1,-1,5), hold on
% plot(CD(:),HD(:),'b.')
% xlim([-.1 0.3])
% ylim([-10 5])
drawnow

end