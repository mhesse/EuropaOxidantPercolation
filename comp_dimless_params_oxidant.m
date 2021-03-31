function [Scales,Param,ic,ConstFuns] = comp_dimless_params_oxidant(Scales,Param,ic)
% author: Marc Hesse
% date: 26 Mar 2019
% Description:
% Takes the characteristic properties of the ice shell and computes the
% characteristic scales and dimensionless govering paramters.

% Additional paramters that are typically not modified
muf = 1e-3;      % viscosity of water [Pa s]
R = 8.314;       % gas constant [J/(mol K)]

Vm = 1.97e-5;    % Molar volume of ice I [m^3/mol]
Dov = 9.1e-4;    % volume diffusion constant [m^2/s]
Qvstar = 59.4e3; % volume diffusion activation energy [J/mol]
Tm = 273.15;     % melting point of water [K]
A = Qvstar/(R*Tm);
Drho = 80;       % density difference between ice and brine [kg/m3]
g = 1.315;       % Europa's surface gravity [m/s^2]

% Permeability model
if strcmp(Param.perm,'vonBargen')
    tau = 1600;      % tortuosity [1]
    Param.exp.n = 2; % Permeability exponent [-]
elseif strcmp(Param.perm,'WarkWatson')
    tau = 200;
    Param.exp.n = 3;
elseif strcmp(Param.perm,'GHP')
    tau = 54;
    Param.exp.n = 3;
else
    error('Unknown permeability model.\n')
end

% Constitutive functions
ConstFuns.k0 = @(d) d.^2/tau;
ConstFuns.k = @(d,phi,n) ConstFuns.k0(d).*phi.^n;
ConstFuns.eta_diff = @(d,T) (R*T.*d.^2)/(42*Vm*Dov).*exp(Qvstar./(R*T));
ConstFuns.eta0 = @(d)            ConstFuns.eta_diff(d,Tm);
ConstFuns.eta  = @(d,T,sw,phi)   ConstFuns.eta0(d).*exp(A*(Tm./T-1)).*exp(-abs(sw)*phi);
ConstFuns.xi   = @(d,T,sw,phi,m) ConstFuns.eta(d,T,sw,phi)./(phi.^m);
ConstFuns.phi_tilde = @(phiD) phiD + (1-Scales.phi_c*phiD).*Param.conc1.Dist/Scales.phi_c;

% Charactristic scales
Scales.k_c   = ConstFuns.k(Scales.d_c,Scales.phi_c,Param.exp.n);   
Scales.eta_c = ConstFuns.eta(Scales.d_c,Scales.T_c,0,Scales.phi_c);
Scales.xi_c  = ConstFuns.xi(Scales.d_c,Scales.T_c,0,Scales.phi_c,Param.exp.m);


Scales.delta   = @(d,T,sw,phi,m,n) sqrt((ConstFuns.xi(d,T,sw,phi,m).*ConstFuns.k(d,phi,n))/muf);
Scales.delta_c = Scales.delta(Scales.d_c,Scales.T_c,0,Scales.phi_c,Param.exp.m,Param.exp.n);
Scales.h_c     = Scales.delta_c;
Scales.p_c     = Drho*g*Scales.delta_c;
Scales.K_c     = Scales.k_c*Drho*g/muf;
Scales.Xi_c    = Scales.xi_c/(Drho*g);
Scales.v_c     = Scales.K_c;
Scales.q_c     = Scales.K_c;
Scales.u_c     = Scales.K_c*Scales.delta_c;
Scales.t_c     = Scales.phi_c*Scales.Xi_c/Scales.delta_c;
Scales.sec2yr  = 60^2*24*365.25;


% Dimensionless paramters
Param.non_dim.Pe = Scales.phi_c;
Param.non_dim.Da = 0;
Param.non_dim.H  = Param.dim.H/Scales.delta_c;
Param.non_dim.L  = Param.dim.L/Scales.delta_c;
Param.non_dim.D  = Param.dim.D/Scales.delta_c;
Param.non_dim.Dhom  = Param.dim.Dhom/Scales.delta_c;
if isfield(Param.dim,'melt')
    if isfield(Param.dim.melt,'thick')
        Param.non_dim.melt.thick = Param.dim.melt.thick/Scales.delta_c;
        Param.dim.melt.above = Param.non_dim.H - Param.non_dim.melt.thick;
    end
end

% Package additional paramters for documentation
Param.dim.muf = muf;       % viscosity of water [Pa s]
Param.dim.R = R;           % gas constant [J/(mol K)]
Param.dim.tau = tau;       % tortuosity [1]
Param.dim.Vm = Vm;         % molar volume of ice I [m^3/mol]
Param.dim.Dov = Dov;       % volume diffusion constant [m^2/s]
Param.dim.Qvstar = Qvstar; % volume diffusion activation energy [J/mol]
Param.dim.Tm = Tm;         % melting point of water [K]
Param.dim.A = A;
Param.dim.Drho = Drho;     % density difference between ice and brine [kg/m3]
Param.dim.g = g;           % Europa's surface gravity [m/s^2]

%% Initial consition
ic.z_sal  = Param.dim.H-ic.d_sal;  % Base of the enriched crust [m]
ic.z_oxy = Param.dim.H-ic.d_oxy;   % Base of oxidant enriched layer [m]
ic.z_per = ic.z_sal-ic.d_per;      % Base of perturbed layer [m]
ic.dD_sal = ic.d_sal/Scales.delta_c;
ic.dD_oxy = ic.d_oxy/Scales.delta_c;
ic.zD_sal = ic.z_sal/Scales.delta_c;
ic.zD_oxy = ic.z_oxy/Scales.delta_c;
ic.zD_per = ic.z_per/Scales.delta_c;

ic.zD_tbl = ic.z_tbl/Scales.delta_c; % depth of the thermal boundary layer