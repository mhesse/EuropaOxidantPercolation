function [dist_CD,dist_HD,dist_cD_tracer] = comp_dist_coeffs(CD,HD,TD,Phi,dist_cD_tracer,reg,rho,hD,Phase)
% authors: Marc Hesse, Jake Jordan
% date: 13 Dec 2019

% Description:
% Uses the binary eutectic phase behavior to compute the non-linear
% distribution coefficients for CD, HD and any tracer

phi_ice = Phi(:,1);
phi_sal = Phi(:,2);
phi_bri = Phi(:,3);

%% Distribution coefficient tracer
dist_cD_tracer = extend_to_single_phase_regions(dist_cD_tracer,reg);

%% Distribution coefficient salt
Xsal = 1;
Xbri = Phase.Xbri(TD);
cD_sal = rho.si*Xsal;
cD_bri = rho.bi*Xbri.*(reg == 4 | reg == 3) + CD.*(reg == 5);
dist_CD = phi_sal./(1-phi_bri).*cD_sal./cD_bri;
dist_CD(reg == 4 & CD == 0) = 0; % deal with pure ice limit in region 4 (ice+bri)
dist_CD = extend_to_single_phase_regions(dist_CD,reg);

%% Distribution coefficient enthalpy
HD_bri = rho.bi*hD.bri(TD).*(reg == 4 | reg == 3)+ HD.*(reg == 5);
HD_ice = hD.ice(TD);
HD_sal = rho.si*hD.sal(TD);
dist_HD = phi_ice./(1-phi_bri).*HD_ice./HD_bri + phi_sal./(1-phi_bri).*HD_sal./HD_bri;
dist_HD = extend_to_single_phase_regions(dist_HD,reg);
end

function [dist_coeff] = extend_to_single_phase_regions(dist_coeff,reg)
dist_coeff(reg == 5) = 0;
dist_coeff(reg == 1 | reg == 2) = 1;
end
