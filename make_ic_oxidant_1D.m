function [CD,HD,TD,Phi,phiD,CD_tracer,cD_tracer,X] = make_ic_oxidant_1D(ic,Scales,Grid,Param,ConstFuns,TD_fun,zDc,HD_of_TX,CD_of_TX,rho,cp,hD,Phase)

% Scale inputs
zD_sal = ic.z_sal/Scales.delta_c;
zD_tbl = ic.z_tbl/Scales.delta_c;

% Mass fraction of salt
X = .5+.5*tanh(5*(zDc-zD_sal));
X = ic.X_con + (ic.X_top-ic.X_con)*X;

% Temperature profile
TD_top = TD_fun(ic.Ttop);
TD_bot = TD_fun(ic.Tbot);
TD = TD_top + (TD_bot-TD_top)/(Grid.ymax-zD_tbl)*(Grid.ymax-zDc);
TD(zDc<zD_tbl) = TD_bot;

% Bulk composition and enthalpy
HD = HD_of_TX(TD,X); HD = HD(:);
CD = CD_of_TX(TD,X); CD = CD(:);
X = X(:);

% Volume fractions
[TD,Phi,regime] = eval_phase_behavior(HD,X,rho,cp,hD,Phase);
phiD = Phi(:,3)/Scales.phi_c;

% Tracer distribution
cD_tracer= zeros(Grid.Ny,1); cD_tracer(zDc>=Param.dim.melt.above) = 1;
CD_tracer = ConstFuns.phi_tilde(phiD).*cD_tracer; 

% Add dimensionless values to structure
% add to structure
ic.zD_sal = zD_sal;
ic.zD_tbl = zD_tbl;
ic.TD_top = TD_top;
ic.TD_bot = TD_bot;
