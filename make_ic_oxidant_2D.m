function [CD,HD,TD,Phi,phiD,CD_tracer,cD_tracer,X] = make_ic_oxidant_2D(ic,Scales,Grid,Param,ConstFuns,TD_fun,ZDc,XDc,HD_of_TX,CD_of_TX,rho,cp,hD,Phase,D,G,I)

% Scale inputs
zD_sal = ic.z_sal/Scales.delta_c; % depth of salt
zD_tbl = ic.z_tbl/Scales.delta_c; % depth of the thermal boundary layer

% Mass fraction of salt
X = .5+.5*tanh(5*(ZDc-zD_sal));       % 1D version
X = ic.X_con + (ic.X_top-ic.X_con)*X;

%% Generate random perturbations in the salt mass fraction
rng_state = ic.seed;
sigma = ic.sigma;
Expert = comp_matern(sigma,ic.corr_length,Grid,D,G,I,rng_state);
X = X + reshape(Expert,Grid.Ny,Grid.Nx);
X(X<0) = 1e-4;
% smooth it a bit
gamma = 1e-2;

for i=1:25
    Im = I - gamma*D*G;
    X = Im\X(:);
    
    X = reshape(X,Grid.Ny,Grid.Nx);

end
%     surf(XDc,ZDc,X)
% zlim([min(X(:)),1e-3])
% caxis([min(X(:)),1e-3])
% drawnow
% min(X(:))

% Temperature profile
TD_top = TD_fun(ic.Ttop);
TD_bot = TD_fun(ic.Tbot);
TD = TD_top + (TD_bot-TD_top)/(Grid.ymax-zD_tbl)*(Grid.ymax-ZDc);
TD(ZDc<zD_tbl) = TD_bot;

% Bulk composition and enthalpy
HD = HD_of_TX(TD,X); HD = HD(:);
CD = CD_of_TX(TD,X); CD = CD(:);
X = X(:);

% Volume fractions
[TD,Phi,regime] = eval_phase_behavior(HD,X,rho,cp,hD,Phase);
phiD = Phi(:,3)/Scales.phi_c;

% Tracer distribution
cD_tracer= zeros(Grid.Ny,1); cD_tracer(ZDc>=Param.dim.melt.above) = 1;
CD_tracer = ConstFuns.phi_tilde(phiD).*cD_tracer; 

% Add dimensionless values to structure
% add to structure
ic.zD_sal = zD_sal;
ic.zD_tbl = zD_tbl;
ic.TD_top = TD_top;
ic.TD_bot = TD_bot;
