function [TD,Phi,field] = eval_phase_behavior(HD,CD,rho,cp,hD,Phase) % repo
% author: Marc Hesse
% date: 19 Jul 2019
% Input:
% HD = N by 1 column vector of dimensionless total enthalpy of the system
% CD = N by 1 column vector of dimensionless total composition of component 2 in system
% rho = structure containing phase densities 
% cp  = structure containing phase heat capacities 
% h   = structure containing the specific enthalpies of the phases
% Phase = structure containing information about phase behavior
%
% Output:
% TD = N by 1 column vector of dimensionless temperatures
% Phi = N by 3 matrix of phase volume fractions [ice,sal,bri]
% field = N by 1 column vector with classification of fieldion in phase diagram [-]
% Description:
% This function computes the phase split for a simple binary eutectic with
% a linear liquidus given the dimensionless total enthalpy, HD, and the
% total composition, CD.


% Separate input into 5 possible fieldions: 
% 1) Pure solid 1:   solid 1
% 2) Sub-solidus:    solid 1 + solid 2
% 3) Eutectic:       solid 1 + solid 2 + brine
% 4) Super-eutectic: solid 1           + brine
% 5) Super-liquidus:                     brine
N = length(HD(:));

% Unpack phase structure for readability
L  = Phase.L;
Xe = Phase.Xe;
T1 = Phase.T1;
TDe = Phase.TDe;
TDl = Phase.TDl;
HDe = Phase.invHC.HDe;
HDl = Phase.invHC.HDl;
CDb = Phase.invHC.CDb;

CD3l = Phase.invHC.CD3l;
Ste = Phase.Ste;

% Identify different fieldions of HX phase diagram
field = zeros(N,1);
field(HD<=1 & CD==0)            = 1;
field(HD<=0 & CD>0)             = 2;
field(HD>0 & HD<=HDe(CD) & CD>0 & CD <= CDb(HD)) = 3;
field((HD>HDe(CD) & HD<=HDl(CD) & CD>0) | (HD>1 & HD<=HDl(0) & CD==0))  = 4;
field(HD>HDl(CD) & CD <= CD3l)               = 5;

% Initialize the output vectors
Phi = zeros(N,3); % [ice,sal,bri]
TD   = zeros(N,1);

%% Regime 1 - Single-phase sub-solius: solid 1 (ice)
% Mass fractions
% Don't need to compute them.

% Volume fractions (no salt and no brine)
Phi(field==1,1) = 1; % ice

% Temperature
TD(field == 1) = HD(field==1); 

%% Regime 2 - Two-phase sub-solidus:    solid 1 + solid 2 (ice + salt)

% Volume fractions (no brine)
Phi(field==2,2) = CD(field==2)/rho.si; % sal
Phi(field==2,1) = 1 - Phi(field==2,2);                  % ice

% Temperature
TD(field == 2) = HD(field==2)./(Phi(field==2,1) + Phi(field==2,2)*rho.si*cp.si); 

%% Regime 3 - Eutectic: solid 1 + solid 2 + melt (ice + salt + brine)
% Volume fraction
Phi(field==3,3) = Ste*HD(field==3)/(rho.bi);                 % bri
Phi(field==3,2) = (CD(field==3)-Ste*Xe*HD(field==3))/rho.si; % sal
Phi(field==3,1) = 1 - Phi(field==3,2) - Phi(field==3,3);     % ice
% Temperature
TD(field == 3) = 0; % Temperature is fixed at the eutectic TDe = 0


%% Regime 4 - Super-eutectic: solid 1 + brine (ice + brine)
% Note 1: The case X=0 (pure ice) is treated separately.
% Note 2: For X > 0, this field requires the solution of a quadratic 
% equation to determine the temperature. The solution to this quadratic 
% is in Phase.field3.T1, which comes from setup_phase_behavior.

% Temperature
if  min(Phase.invHC.field4.disc(CD(field==4 & CD(:)~=0),HD(field==4 & CD(:)~=0)) < 0) error('Negative discriminant!'); end

TD(field == 4 & CD(:)~=0) = Phase.invHC.field4.T2(CD(field==4 & CD(:)~=0),HD(field==4 & CD(:)~=0));
TD(field == 4 & CD(:)==0) = 1; % single phase limit of fieldion 4

% Volume fractions (no salt)
Phi(field==4 & CD(:)~=0,3) =  CD(field == 4 & CD(:)~=0)./(rho.bi*Xe*(1-TD(field == 4 & CD(:)~=0))); % brine
Phi(field==4 & CD(:)==0,3) = (HD(field == 4 & CD(:)==0)-1)/(rho.bi/Ste+rho.bi*cp.bi-1); % single phase limit
Phi(field==4,1) = 1 - Phi(field==4,3); % ice


%% Regime 5 - Super-liquidus: brine
% Temperature
TD(field==5) = TDe+(HD(field==5)-rho.bi/Ste)./(rho.bi*cp.bi);
% Volume fractions
Phi(field==5,3) = 1;

%% Unidentified 
N_unid = sum(field == 0);
if N_unid > 0
    error('%d unidentified points in HC-space!\n',N_unid)
    TD(field == 0) = nan;
    Phi(field == 0,:) = nan;
end
Phi(Phi<0)=0;