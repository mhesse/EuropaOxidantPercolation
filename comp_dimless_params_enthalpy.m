function [Scales,Param,TD_fun,HD_fun] = comp_dimless_params_enthalpy(Scales,Param,Phase,rho,cp,kappa)
% TeD = Phase.TeD;
TD_fun = @(T) (T-Phase.Te)/(Phase.T1-Phase.Te);
HD_fun = @(Phi,T) (1-Phi(:,3)).*TDic+phi_bri*rho.bi.*(1/Phase.Ste+cp.bi*TDic);
Scales.kappa_c = kappa.ice;%(Scales.T_c);
Scales.D_c = Scales.kappa_c/(rho.ice*cp.ice);
Scales.t_T = Scales.delta_c^2/Scales.D_c;

% Dimensionless parameters
% Conductive Peclet number
% Can also be interpreted as Pe = ve_c*x_c/D_c, where ve_c is the
% characteristic effective velocity in the enthalpy equation.
% See notes from 25 Dec 2019 for details.
Param.non_dim.PeT = Scales.delta_c^2/(Scales.D_c*Scales.t_c); 

