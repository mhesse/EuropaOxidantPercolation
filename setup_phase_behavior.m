function [HD,CD,phi,f,Phase] = setup_phase_behavior(rho,cp,hD,Phase) % repo
% author: Marc Hesse
% date: 1 May 2019
% Input:
% rho = structure containing phase densities (both dimensional and dimensionless)
% cp  = structure containing phase heat capacities (both dimensional and dimensionless)
% hD  = structure containing the dimensionless specific enthalpies of the phases
% Phase = structure containing information about phase behavior
%
% Output:
% HD = anonymous function for dimensionless total enthalpy in terms of X and TD
% CD = anonymous function for dimensionless total composition in terms of X and TD
% phi = structure containing volume fraction as function of dimensionless 
%       total enthalpy and composition, HD and CD.
% f   = structure containing mass fraction as function of dimensionless 
%       total enthalpy and composition, HD and CD.
% Phase = updated structure with quadratic solution for supra-eutectic
% region
%
% Description:
% This function does several thing. 
% 1) First, it provides functions that compute other quantities 
%    (HD,CD, phi's, f's) given the native variables, X and TD, of the phase 
%    diagram. 
% 2) Second, it provides the boundaries of the different fields in the 
%    phase diagram in different variables. 
% 3) Third, it gives the solution of the quadratic required in the 
%    supra-eutectic two-phase field in the HX and HC spaces.

% Unpack phase structure for ease of reading
L  = Phase.L;
Xe = Phase.Xe;
T1 = Phase.T1;
Te = Phase.Te;
Ste = Phase.Ste;
DT = T1-Te;

% Define liquidus and phase composition
TDe = 0;
TD1 = 1;
TDl = @(X) 1-X/Xe;
Xbri = @(TD) Xe*(1-TD);


% Mass fractions
f.sal = @(TD,X) X.*(TD<TDe);
f.bri = @(TD,X) (X./Xbri(TD)).*(TD<=TDl(X) & TD>=TDe) + (TD>TDl(X));
f.ice = @(TD,X) (1-f.sal(TD,X)).*(TD<TDe) + ...
                (1-f.bri(TD,X)).*(TD<=TDl(X) & TD>=TDe);

% Volume fractions
phi.sal = @(TD,X) f.sal(TD,X)./(rho.si + (1-rho.si)*f.sal(TD,X));
phi.bri = @(TD,X) f.bri(TD,X)./(rho.bi + (1-rho.bi)*f.bri(TD,X));
phi.ice = @(TD,X) (1-phi.sal(TD,X)).*(TD<TDe) + ...
                  (1-phi.bri(TD,X)).*(TD<=TDl(X) & TD>=TDe);

% dimensionless CD as function of TD and X
CD = @(TD,X) phi.sal(TD,X).*rho.si+...
             phi.bri(TD,X).*rho.bi.*Xbri(TD);

% dimensionless HD as function of TD and X
HD = @(TD,X) phi.ice(TD,X).*hD.ice(TD)+...
             phi.sal(TD,X).*rho.si.*hD.sal(TD)+...
             phi.bri(TD,X).*rho.bi.*hD.bri(TD);

%% 2a) Boundaries of HX - phase fields 
% Boundaries of the regions in HX-phase diagram
Phase.invHX.HDe = @(X) rho.bi/Ste*X./(rho.bi*Xe+(1-rho.bi)*X);
Phase.invHX.HDl = @(X) rho.bi*(1/Ste+cp.bi*(1-X/Xe));

% Corners of phase fields in HC-diagram
% 1
Phase.invHX.X1  = 0;
Phase.invHX.HD1 = 0;

% 2s
Phase.invHX.X2s  = 0;
Phase.invHX.HD2s = 1;
% 2l
Phase.invHX.X2l  = 0;
Phase.invHX.HD2l = rho.bi*(1/Phase.Ste + cp.bi);

% 3s
Phase.invHX.X3s  = Xe;
Phase.invHX.HD3s = 0;
% 3l
Phase.invHX.X3l  = Xe;
Phase.invHX.HD3l = rho.bi/Phase.Ste;

% Test consistency between boundaries and corners
if Phase.invHX.HD2l - Phase.invHX.HDl(0) >1e-13; error('Left endpoint of HDl does not match!\n'); end
if Phase.invHX.HD3l - Phase.invHX.HDl(Xe) > 1e-13; error('Right endpoint of HDl does not match!\n'); end

%% 2b) Boundaries of HC - phase fields 
% Boundaries of the regions in HC-phase diagram
nu = (1-rho.si)*Phase.Xe+rho.si;
Phase.invHC.HDe = @(CD) CD/(Xe*Ste);
Phase.invHC.HDl = @(CD) rho.bi*(1/Ste+cp.bi*(1-CD/(rho.bi*Xe)));
Phase.invHC.HDb = @(CD) rho.bi*rho.si/((rho.bi*nu-rho.si)*Ste)*(nu*CD/(rho.si*Xe)-1);
Phase.invHC.CDb = @(HD) rho.si*Xe/nu*(1-Ste*HD/rho.bi)+Ste*Xe*HD;
Phase.invHC.nu = nu;

% Corners of phase fields in HC-diagram
% 1
Phase.invHC.CD1 = 0;
Phase.invHC.HD1 = 0; 

% 2s
Phase.invHC.CD2s = 0;
Phase.invHC.HD2s = 1;
% 2l
Phase.invHC.CD2l = 0;
Phase.invHC.HD2l = Phase.invHX.HD2l;

% 3s
Phase.invHC.CD3s = rho.si*Xe/((1-rho.si)*Xe+rho.si);
Phase.invHC.HD3s = 0;
% 3l
Phase.invHC.CD3l = rho.bi*Xe;
Phase.invHC.HD3l = Phase.invHX.HD3l;
      
%% 3) Solution to quadratic in supra-eutectic region
% Solution in HX-phase diagram
a = @(X)  rho.bi*Xe./X+(1-rho.bi);
b = @(X) -rho.bi*Xe./X;

alpha = @(X) b(X);
beta = @(X,HD) a(X)+rho.bi*cp.bi-1-b(X).*HD;
gamma = @(X,HD) rho.bi/Ste-a(X).*HD;

% Roots
disc_HX = @(X,HD) beta(X,HD).^2 - 4*alpha(X).*gamma(X,HD);

Phase.invHX.field4.T1 = @(X,HD) (-beta(X,HD) + sqrt(disc_HX(X,HD)))./(2*alpha(X));
Phase.invHX.field4.T2 = @(X,HD) (-beta(X,HD) - sqrt(disc_HX(X,HD)))./(2*alpha(X));

Phase.invHX.field4.disc = disc_HX;
Phase.invHX.field4.alpha = alpha;
Phase.invHX.field4.beta = beta;
Phase.invHX.field4.gamma = gamma;

% Solution HC-phase diagram
alpha_HC = rho.bi*Xe;
beta_HC  = @(CD,HD) (1-rho.bi*cp.bi)*CD - (1+HD)*rho.bi*Xe;
gamma_HC = @(CD,HD) rho.bi*Xe*HD - rho.bi/Ste*CD;
disc_HC = @(CD,HD) beta_HC(CD,HD).^2 - 4*alpha_HC.*gamma_HC(CD,HD);

Phase.invHC.field4.T1 = @(CD,HD) (-beta_HC(CD,HD) + sqrt(disc_HC(CD,HD)))./(2*alpha_HC);
Phase.invHC.field4.T2 = @(CD,HD) (-beta_HC(CD,HD) - sqrt(disc_HC(CD,HD)))./(2*alpha_HC);

Phase.invHC.field4.disc = disc_HC;
Phase.invHC.field4.alpha = alpha_HC;
Phase.invHC.field4.beta = beta_HC;
Phase.invHC.field4.gamma = gamma_HC;

%% 4) Clean-up: Add other functions to Phase structure

Phase.TDl = TDl;
Phase.Xbri = Xbri;
Phase.TDe = TDe;
Phase.TD1 = TD1;
Phase.DT = DT;

