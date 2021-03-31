function [CD,HD,TD,Phi,phiD,CD_tracer,cD_tracer,X] = make_ic_oxidant_2D(ic,Scales,Grid,Param,ConstFuns,TD_fun,ZDc,XDc,HD_of_TX,CD_of_TX,rho,cp,hD,Phase,D,G,I)


% Mass fraction of salt
X = .5+.5*tanh(5*(ZDc-ic.zD_sal));       % 1D version
X = ic.X_con + (ic.X_top-ic.X_con)*X;

%% Generate random perturbations in the salt mass fraction
% X = marc_fancy_perturb(X,ic,Grid,D,G,I);
ind = ZDc(:) < ic.zD_sal & ZDc(:) > ic.zD_per;
ampl = ic.X_con/10;
RAMP = ampl*(ZDc(ind)-ic.zD_per)/(ic.zD_sal-ic.zD_per);
X(ind) = ic.X_con + RAMP.*cos(ic.n_pert*pi*XDc(ind)/Grid.Lx);
% surf(XDc,ZDc,log10(reshape(X,Grid.Ny,Grid.Nx)))
% ylim([100 ic.zD_sal])

% Temperature profile
TD_top = TD_fun(ic.Ttop);
TD_bot = TD_fun(ic.Tbot);
TD = TD_top + (TD_bot-TD_top)/(Grid.ymax-ic.zD_tbl)*(Grid.ymax-ZDc);
TD(ZDc<ic.zD_tbl) = TD_bot;

% Bulk composition and enthalpy
HD = HD_of_TX(TD,X); HD = HD(:);
CD = CD_of_TX(TD,X); CD = CD(:);
X = X(:);

% Volume fractions
[TD,Phi,regime] = eval_phase_behavior(HD,X,rho,cp,hD,Phase);
phiD = Phi(:,3)/Scales.phi_c;

% Tracer distribution
cD_tracer= zeros(Grid.Ny,1); 
if strcmp(ic.oxy_mix,'no')
    cD_tracer(ZDc>=ic.zD_oxy) = 1;
elseif strcmp(ic.oxy_mix,'yes')
    cD_tracer(ZDc>=ic.zD_sal) = ic.dD_oxy/ic.dD_sal;
    fprintf('Oxidant is mixed.\n')
else
    error('Unknown option in ic.oxy_mix.')
end
CD_tracer = ConstFuns.phi_tilde(phiD).*cD_tracer; 

% Add dimensionless values to structure
% add to structure
ic.TD_top = TD_top;
ic.TD_bot = TD_bot;
end

function [X] = marc_fancy_perturb(X,ic,Grid,D,G,I)
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
end