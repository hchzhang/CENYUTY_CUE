%clear all;
%close all;
function kk_out = kk_matrix(ts,sp,tem,swc)
% This function was used to calculate the decomposition rate of each carbon
% pool based on the turnover time, temperature control, and moisture
% control

global np   % fix seting
global clay lignin_struc_cmatrix  % for sensitivity test
global tau_metabolic tau_struct tau_active tau_slow tau_passive flux_tot_coeff
global litter_struct_coef

kk_matrix(1:np,1:np) = 0.0 ;

iabove =1;
ibelow=2;
imetabolic=1;
istructural= 2;

litter_tau(imetabolic) = tau_metabolic  ;
litter_tau(istructural) = tau_struct ;

soc_tau = [tau_active,tau_slow,tau_passive]; % residence time in active, slow, passive pools (years)

frozen_respiration_func = 1 ;
if (ts>sp)
    tsurf_in =tem+273.15; % temperature of aboveground litter (land surface )
    tsoil_decomp = tem+273.15; % temperature of belowground litter (deep soil)
else
    tsurf_in =tem+273.15; % temperature of aboveground litter (land surface )
    tsoil_decomp = tem+273.15; % temperature of belowground litter (deep soil)
end
% Incubated litter
% if (ts>sp+15+77)
%     tsurf_in =295.2; % temperature of aboveground litter (land surface )
%     tsoil_decomp = 295.2; % temperature of belowground litter (deep soil)
% elseif (ts>sp+77)
%     tsurf_in =295.2; % temperature of aboveground litter (land surface )
%     tsoil_decomp = 295.2; % temperature of belowground litter (deep soil)
% else
%     tsurf_in =298.2; % temperature of aboveground litter (land surface )
%     tsoil_decomp = 298.2; % temperature of belowground litter (deep soil)
% end

litterhum = swc; % humidity of aboveground litter 
soilhum_decomp = swc; % humidity of belowground litter 

[control_temp_som(iabove),control_temp_litt(iabove)] = control_temp_func (tsurf_in, frozen_respiration_func,tem);
[control_temp_som(ibelow),control_temp_litt(ibelow)] = control_temp_func (tsoil_decomp, frozen_respiration_func,tem);
control_moist(iabove) = control_moist_func(litterhum);
control_moist(ibelow) = control_moist_func(soilhum_decomp);
%============================================================================

kk_matrix(1,1) = 1.0/litter_tau(imetabolic) * control_temp_litt(iabove)...
    * control_moist(iabove);
kk_matrix(2,2) = 1.0/litter_tau(imetabolic) * control_temp_litt(ibelow)...
    * control_moist(ibelow);
kk_matrix(3,3) = 1.0/litter_tau(istructural) * control_temp_litt(iabove)...
    * control_moist(iabove)* exp( -litter_struct_coef * lignin_struc_cmatrix(1) );
kk_matrix(4,4) = 1.0/litter_tau(istructural) * control_temp_litt(ibelow) ...
    * control_moist(ibelow)* exp( -litter_struct_coef * lignin_struc_cmatrix(2) );

kk_matrix(5,5)= 1.0/soc_tau(1) * control_moist(ibelow) ...
    * control_temp_som(ibelow) * ( 1. - flux_tot_coeff(3) * clay);
kk_matrix(6,6)= 1.0/soc_tau(2) * control_moist(ibelow)...
    * control_temp_som(ibelow);
kk_matrix(7,7)= 1/soc_tau(3)*control_moist(ibelow)...
    * control_temp_som(ibelow);

kk_out = kk_matrix;

end
