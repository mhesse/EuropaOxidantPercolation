function [rho,cp,kappa,hD,HD,Phase] = physical_properties_eutectic(Phase,Scales) % repo
% author: Marc A. Hesse
% date: 1,9 May 2019

% DESCRIPTION:
% Defines the physical properties of the phases 

% INPUT:
% Te = eutectic temperature [K]
% T1 = melting temperature of solid 1 [K]
% T2 = melting temperature of solid 2 [K]
% L = latent heat of water [J/kg]
% phi = structure containing the initial volume fractions of phases [-]
%       phi.ice0, phi.hyd0, phi.sal0, phi.sil0

% OUTPUT:
% phi = structure containing the function handles for the volume fractions 
%       of the phases as function of temperature [-]
% rho = structure containing the constant densities of the phases [kg/m^3]
% cp = structure containing the function handles for the specific heat 
%      capacities of the phases as function of temperature [J/(kg K)]
% kappa = structure containing the function handles for the thermal 
%      conductivities of the phases as function of temperature [W/(m K)]
% hD = structure containing the function handles for the dimensionless 
%      specific enthalpy of the phases as function of temperature [-]
% HD = function handle for dimensionless bulk entalpy of system [-]

% EXAMPLE CALL
% >> Te = 245;
% >> T1 = 273;
% >> T2 = 400;
% >> L = 3.34e5;
% >> phi0.ice = 0.3;
% >> phi0.hyd = 0.3;
% >> phi0.sal = 0.2;
% >> [phi,rho,cp,kappa,h,chi] = physical_properties_eutectic(Ts,Tl,L,phi0);

% Unpack phase structure
L  = Phase.L;
Xe = Phase.Xe;
T1 = Phase.T1;
Te = Phase.Te;
DT = T1-Te; Phase.DT = DT;

%% Densities of the phases 
% % dimensional [kg/m^3]
rho.ice = 917;  % ice
rho.sal = 1466; % Mirabillite
rho.bri = 1010; % brine
% Busynesk approximation!
% rho.ice = 1e3;  % ice
% rho.sal = 1e3; % Mirabillite
% rho.bri = 1e3; % brine
if rho.sal == rho.ice && rho.ice == rho.bri
disp('Densities adjusted for Boussinesq approximation!')
end
% dimensionless density ratios [1]
rho.bi = rho.bri/rho.ice;
rho.si = rho.sal/rho.ice;

%% Heat capacity 
% dimensional [J/kg/K]
cp.ice = 2000;
cp.sal = 920;
cp.bri = 4200;
% dimensionless heat capacity ratios [1]
cp.bi  = cp.bri/cp.ice;
cp.si  = cp.sal/cp.ice;
Ste = cp.ice*DT/L;
Phase.Ste = Ste;

%% Thermal conductivity [W/m/K]
kappa.ice = @(T) 0.4685 + 488.12./T;
kappa.sal = @(T) 0.6  + T*0;
kappa.bri = @(T) 0.56 + T*0;
% contant properties
kappa.ice = kappa.ice(262);
kappa.sal = kappa.sal(262);
kappa.bri = kappa.bri(262);

kappa.ii = kappa.ice/kappa.ice;
kappa.si = kappa.sal/kappa.ice;
kappa.bi = kappa.bri/kappa.ice;

kappa.sysD = @(TD,Phi) Phi(:,1).*kappa.ii+...
                       Phi(:,2).*kappa.si+...
                       Phi(:,3).*kappa.bi;
             

%% Dimensionless Specific enthalpy of the phases
hD.ice = @(TD) TD;
hD.sal = @(TD) cp.si*TD;
hD.bri = @(TD) hD.ice(0) + 1/Ste + cp.bi*TD;

%% Bulk Enthalpy of system
HD = @(TD,Phi) Phi(:,1).*hD.ice(TD)+...
               Phi(:,2).*rho.si.*hD.sal(TD)+...
               Phi(:,3).*rho.bi.*hD.bri(TD);
