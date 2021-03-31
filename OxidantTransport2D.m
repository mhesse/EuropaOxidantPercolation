% OxidantTransport.m
% author: Marc Hesse, Jake Jordan
% date: 23 Mar 2020
% Description:
% This code solves the Europa oxidant transport problem (high pororsity and 
% tracer at top) with heat conduction and phase behavior in 2D. The version 
% before this was CH_OxidantTransport2D_variable_dt. This version includes 
% the rejection of the oxidant from the solid region.
close all, clear, clc

%% Define physical parameters

% Domain Geometry
% Param.refine = 4;              % Obsolete
Param.dim.H = 25e3;              % Ice shell thickness  [m]
Param.dim.L = 3e3;               % Width of the domain  [m]
Param.dim.D = 5e3;               % Depth of ocean below [m]
Param.dim.Dhom = 2e3;            % Depth below we homogenize composition of ocean [m]
Param.exp.m = 1;                 % compaction viscosity exponent [-]
Param.conc1.Dist = 0e-1;         % Dist = 0 -> all in fluid
Param.phi_lim = 1e-5;            % phi threshold in solid regions
% Param.perm = 'vonBargen';
% Param.perm = 'WarkWatson';
Param.perm = 'GHP';
Param.limiter = 'mc';

% Initial condition
ic.n_pert = 1;            % wave number of perturbation  
ic.d_sal = 3e3;           % depth of salt enrichment [m]
ic.d_oxy = 10e2;          % depth of oxidant enrichment [m]
ic.d_per = 1e3;           % depth of perturbation [m]
ic.z_tbl  = 25e3;         % Base of the thermal boundary layer  [m]
ic.oxy_mix = 'yes';       % oxidant is mixed over whole near surface melt region [m]
ic.X_top = 0.01/0.263;    % Salt (HH) mass fraction in enriched crust [1]
ic.X_con = 1e-4/0.263;    % Salt mass fraction in convecting shell [1]
ic.Ttop = 100;            % Mean surface temperature [K]
ic.Tbot = 262;            % Basal temperature [K]
ic.seed = 4;              % seed for random perturbations
ic.corr_length = 5;       % correlation length of initial perturbations
ic.sigma = ic.X_con*5e-5; % standard deviation of perturbation

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

[Scales,Param,ic,ConstFuns] = comp_dimless_params_oxidant(Scales,Param,ic);
[Scales,Param,TD_fun] = comp_dimless_params_enthalpy(Scales,Param,Phase,rho,cp,kappa);

%% Dimensionless parameters
Param.non_dim.tmax = 25;
Param.non_dim.Pe_oxy = 1e-6;
TD_top = TD_fun(ic.Ttop);
TD_bot = TD_fun(ic.Tbot);

%% Define Grid and build operators
% Grid.periodic = 'x-dir';
Grid.geom = 'cylindrical_rz';
Grid.xmin = 0;  Grid.xmax = Param.non_dim.L; Grid.Nx = 60; %Param.refine*ceil(Param.non_dim.L);
Grid.ymin = -Param.non_dim.D;  Grid.ymax = Param.non_dim.H; Grid.Ny = 600; %Param.refine*ceil(Param.non_dim.H+Param.non_dim.D);
Nt = Param.non_dim.tmax*2e3; Param.num.Nt = Nt;
Param.plot_interval = 100;

theta = 1; Param.num.theta = theta;% forward  Euler (explicit) - Note: going implicit screws up conservation

Grid    = build_grid(Grid);
[D,G,I] = build_ops(Grid);
[XDc,ZDc] = meshgrid(Grid.xc,Grid.yc);
[XDx,ZDx] = meshgrid(Grid.xf,Grid.yc);
[XDz,ZDz] = meshgrid(Grid.xc,Grid.yf);
dtD     = Param.non_dim.tmax/Nt; Param.num.dtD = dtD; 
% zDc = Grid.yc; zDf = Grid.yf;
dof_hom = Grid.dof(ZDc(:)<=-Param.non_dim.Dhom);

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
phiD(abs(phiD)<1e-12) = 1e-5; % just for initial perturbations
[headD,pD,qD] = solve_Helmholtz(D,G,I,phiD,Param.exp.n,Param.exp.m,Grid,B_h,N_h,fn_h,Param_h,ZDc);
[uD,vD]    = solve_Poisson(D,G,I,phiD,Param.exp.m,pD,Grid,B_u,N_u,fn_u,Param_u);
h0 = figure('units','normalized','outerposition',[0 0 .5 1]);
plot_solution_CH_2D(Grid,phiD,cD_tracer,CD,HD,TD,CD_tracer,hD,uD,qD,vD,pD,ZDc,XDc,Scales,Phase,0)

save(['solution_0.mat'],'tD','i','k','phiD','X','CD','HD','headD','uD','qD','vD','pD','CD_tracer','cD_tracer','TD','Phi','Grid','Scales','Phase','ic','XDc','ZDc','Param')

% Time evolution variables
phi_total = zeros(Nt+1,1); phi_total(1) = sum(phiD*Grid.dy);
Ox_total   = zeros(Nt+1,1);   Ox_total(1) = sum(CD_tracer*Grid.dy);
Ox_ocean = zeros(Nt+1,1);
time = zeros(Nt+1,1);
Pe = Param.non_dim.Pe;
dtD_ref = dtD; 
i = 0;
% h3 = figure;
% h1 = figure('units','normalized','outerposition',[.5 .5 .25 .25]);
% semilogy([0,Param.non_dim.tmax],dtD*[1 1],'b-'), hold on
% xlabel 'tD', ylabel 'dtD'

while tD < Param.non_dim.tmax
    i = i+1;
    %% find time step
    tD = tD + dtD; %dtD*i;
    
    %% Phase behavior and distribution coefficients
    HD(dof_hom) = (HD(dof_hom)'*Grid.V(dof_hom))/sum(Grid.V(dof_hom)); 
    CD(dof_hom) = (CD(dof_hom)'*Grid.V(dof_hom))/sum(Grid.V(dof_hom)); 
    CD_tracer(dof_hom) = (CD_tracer(dof_hom)'*Grid.V(dof_hom))/sum(Grid.V(dof_hom));
    [TD,Phi,reg] = eval_phase_behavior(HD,CD,rho,cp,hD,Phase);
    [dist_CD,dist_HD,dist_tracer] = comp_dist_coeffs(CD,HD,TD,Phi,dist_tracer,reg,rho,hD,Phase);
    
    % Temperature/Enthalpy BC's
    TD(Param_H.dof_dir) = [TD_bot*ones(Grid.Nx,1);TD_top*ones(Grid.Nx,1)];
    HD(Param_H.dof_dir) = [HD_from_TPhi(TD_bot*ones(Grid.Nx,1),Phi(Grid.dof_ymin,:));...
                           HD_from_TPhi(TD_top*ones(Grid.Nx,1),Phi(Grid.dof_ymax,:))];
    
    %% Reject oxidant from solid region
    dof_solid = Grid.dof(reg == 1 | reg==2);
    if ~isempty(dof_solid)
        [dof_solid_f_bnd,dof_solid_f] = find_faces(dof_solid,D,Grid);
        dof_molten = Grid.dof; dof_molten(dof_solid) = [];
        dof_solid_f_bnd = dof_solid_f_bnd(dof_solid_f_bnd>Grid.Nfx); % Exclude x-faces (due to cyl-coords)
        dof_solid_f_bnd = setdiff(dof_solid_f_bnd,[Grid.dof_f_xmin;Grid.dof_f_xmax;Grid.dof_f_ymax]);
        [X_f_bnd,Y_f_bnd] = comp_face_coords(dof_solid_f_bnd,Grid); z_front = mean(Y_f_bnd(:));
        [dof_bnd_in,dof_bnd_out] = find_bnd_cells(dof_solid,dof_molten,dof_solid_f_bnd,D,Grid);
        if length(dof_bnd_in) == length(dof_bnd_out)
            CD_tracer(dof_bnd_out) = CD_tracer(dof_bnd_out) + CD_tracer(dof_bnd_in);
            CD_tracer(dof_bnd_in) = 0; cD_tracer(dof_bnd_in) = 0;
        else
            error('Solidification boundary is not straight.\n')
        end
        G_oxy = G; G_oxy(dof_solid_f_bnd,:) = 0;
%         figure(h3), clf
%         plot(XDc(dof_solid),ZDc(dof_solid),'ro'), hold on
%         plot(X_f_bnd,Y_f_bnd,'b-','linewidth',1)
%         plot(XDc(dof_bnd_out),ZDc(dof_bnd_out),'go')
%         plot(XDc(dof_bnd_in),ZDc(dof_bnd_in),'g+')
    else
        G_oxy = G;
    end
    
%     top = find(reg,1,'last'); % top of melting region
%     CD_tracer(top) = CD_tracer(top)+sum(CD_tracer(Phi(:,3)==0));
%     CD_tracer(Phi(:,3)==0) = 0;                   
                       
    % Compute total mass of fluid and tracer
    phiD = Phi(:,3)/Scales.phi_c;
    phi_total(i+1) = sum(phiD*Grid.dx*Grid.dy); C_total(i+1) =  sum(CD*Grid.dx*Grid.dy); time(i+1) = tD;
    
    % Compute tracer in ocean
    Ox_ocean(i+1) =  sum(CD_tracer(ZDc(:)<0).*Grid.V(ZDc(:)<0));

    %% Flow calculations
    phiD(abs(phiD)<1e-12) = Param.phi_lim; 
    [headD,pD,qD] = solve_Helmholtz(D,G,I,phiD,Param.exp.n,Param.exp.m,Grid,B_h,N_h,fn_h,Param_h,ZDc);
    [uD,vD]    = solve_Poisson(D,G,I,phiD,Param.exp.m,pD,Grid,B_u,N_u,fn_u,Param_u);
    
    %% Transport calculations
    % Update porosity
    Aq = flux_upwind(qD,Grid);
    Av = flux_upwind(vD,Grid);
%     Aq_tracer = build_adv_op(qD,CD_tracer,dtD,G,Grid,Param_c,Param.limiter);
%     Av_tracer = build_adv_op(vD,CD_tracer,dtD,G,Grid,Param_c,Param.limiter);
%     
    Ae_CD = comp_effective_flux(phiD,Aq,Av,dist_CD,Grid,Scales);
    Ae_HD = comp_effective_flux(phiD,Aq,Av,dist_HD,Grid,Scales);
%     Ae_cD = comp_effective_flux(phiD,Aq_tracer,Av_tracer,dist_tracer,Grid,Scales);
    Ae_cD = comp_effective_flux(phiD,Aq,Av,dist_tracer,Grid,Scales);
    ve_cD = Ae_cD*ones(Grid.N,1);
    Ae_cD_lim = build_eff_adv_op(Ae_cD,ve_cD,CD_tracer,dtD,G,Grid,Param_c,Param.limiter);
    
    [dtDmax] = comp_time_step(phiD,qD,vD,XDc,ZDc,XDx,ZDx,XDz,ZDz,Grid,ve_cD);
    Ex_CD = I - dtD*D*Ae_CD;
    Ex_cD = I - dtD*D*(Ae_cD_lim-Param.non_dim.Pe_oxy*G_oxy);
    CD = solve_lbvp(I,Ex_CD*CD,B_C,Param_C.g,N_C);
    CD_tracer = solve_lbvp(I,Ex_cD*CD_tracer,B_c,Param_c.g,N_c);
    cD_tracer = CD_tracer./ConstFuns.phi_tilde(phiD);
    if ~isempty(dof_solid) % necessary because due to the low porosity fudge there is flow across the solidification bnd
        CD_tracer(dof_bnd_in) = 0; cD_tracer(dof_bnd_in) = 0;
    end
    % Direct update of HD in mixed HD-TD formulation
    el_kap = reshape(kappa.sysD(TD,Phi),Grid.Ny,Grid.Nx);
    Kappa = comp_mean(el_kap,-1,1,Grid);
    HD = HD - dtD*D*(Ae_HD*HD - 1/Param.non_dim.PeT*Kappa*G*TD);    
    
    if mod(i,Param.plot_interval) == 0
        k = k+1;
        fprintf('total: i = %d: tD= %3.2f;\n',i,tD);
        figure(h0)
        plot_solution_CH_2D(Grid,phiD,cD_tracer,CD,HD,TD,CD_tracer,hD,uD,qD,vD,pD,ZDc,XDc,Scales,Phase,z_front)
        save(['solution_',num2str(k),'.mat'],'tD','i','k','phiD','CD','HD','hD','uD','qD','vD','pD','CD_tracer','cD_tracer','TD','Phi')
        
%         figure(h1)
%             plot(tD,dtDmax.q,'r.',tD,dtDmax.v,'b.',tD,dtDmax.ve,'g.',tD,dtDmax.tot,'ko'), hold on
        dtDmax;
        drawnow
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

fprintf('\nOxidantTransport.m has finished.\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_solution_CH_2D(Grid,phiD,cD,CD,HD,TD,CD_tracer,hD,uD,qD,vD,pD,ZDc,XDc,Scales,Phase,z_front)
T = Phase.Te + (Phase.T1-Phase.Te)*TD;

subplot 131
surf(XDc,ZDc,Scales.phi_c*reshape(phiD,Grid.Ny,Grid.Nx)), view(2), shading interp, colorbar
xlabel 'x_D'
ylabel 'z_D'
title '\phi_D'
axis equal tight

subplot 132
surf(XDc,ZDc,reshape(cD,Grid.Ny,Grid.Nx)), view(2), shading interp, colorbar
xlabel 'x_D'
ylabel 'z_D'
title 'c_D oxidant'
axis equal tight
caxis([0 3])
% plot(cD(Grid.dof_xmin),Grid.yc,'.-'), hold on
% plot([0 max(cD(Grid.dof_xmin))],z_front*[1 1],'k-'), hold off
% ylim([Grid.ymin,Grid.ymax])
subplot 133
surf(XDc,ZDc,reshape(CD_tracer,Grid.Ny,Grid.Nx)), view(2), shading interp, colorbar
xlabel 'x_D'
ylabel 'z_D'
title 'C_D = \phi_D*c_D'
axis equal tight
caxis([0 300])


% plot_HC_diagram(Phase,1,-1,5), hold on
% plot(CD(:),HD(:),'b.')
% xlim([-.1 0.3])
% ylim([-10 5])
drawnow

end