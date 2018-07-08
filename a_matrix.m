function [a_out,nflux_ma,nmi_ma] = a_matrix(ts,sp,matrix_c,ini_litt,mineral_n,n_cue)
% This function is used to calculate the transfer ratio between each of the
% two carbon/litter pools.

global np emax  % fixed seting
global clay frac_soil_metab_aa frac_soil_metab_ab frac_soil_struct_aa frac_soil_struct_ab ...
    frac_soil_struct_sa frac_soil_struct_sb lignin_struc_cmatrix
global frac_passive_active frac_active_slow frac_passive_slow frac_active_passive frac_slow_passive
global cn_r

a_matrix = zeros(np,np);
nflux_ma = zeros(np,np);
nmi_ma = zeros(np,np);

% C decomposition
for j=1:np
    a_matrix(j,j) = -1.0;
end
% C flux between two C/litter pools
% CUE-mineral
%%%%  Formula default  %%%%%
% a_matrix(5,1) = x(1) ;    % above metabolic to active soil
% a_matrix(5,2) = x(1) ;    % below metabolic to active soil
% a_matrix(5,3) = x(2) * (1- lignin_struc_cmatrix(1)) ;  % above structural to active soil
% a_matrix(5,4) = x(2) * (1 - lignin_struc_cmatrix(2)) ; % below structural to active soil
% a_matrix(6,3) = x(2) * lignin_struc_cmatrix(1); % above structural to slow soil
% a_matrix(6,4) = x(2) * lignin_struc_cmatrix(2); % below structural to slow soil
%%%%  Formula 1  %%%%%
% a_matrix(5,1) = emax * min(1,power(max(1,cn_r(1)/cn_r(5)),x(1)*(mineral_n-x(2)))) ;  % above metabolic to active soil
% a_matrix(5,2) = emax * min(1,power(max(1,cn_r(2)/cn_r(5)),x(1)*(mineral_n-x(2)))) ;  % below metabolic to active soil
% % a_matrix(5,3) = emax * min(1,power(max(1,cn_r(3)/cn_r(5)),x(1)*(mineral_n-x(2))))*(1- lignin_struc_cmatrix(1)) ;  % above structural to active soil
% % a_matrix(5,4) = emax * min(1,power(max(1,cn_r(4)/cn_r(5)),x(1)*(mineral_n-x(2))))*(1 - lignin_struc_cmatrix(2)) ; % below structural to active soil
% % a_matrix(6,3) = emax * min(1,power(max(1,cn_r(3)/cn_r(6)),x(1)*(mineral_n-x(2))))*lignin_struc_cmatrix(1); % above structural to slow soil
% % a_matrix(6,4) = emax * min(1,power(max(1,cn_r(4)/cn_r(6)),x(1)*(mineral_n-x(2))))*lignin_struc_cmatrix(2); % below structural to slow soil
% %%%%  Formula 2  %%%%%
% thres1=x(2)*cn_r(1)/cn_r(5)*ini_litt;
% a_matrix(5,1) = emax * min(1,power(max(1,cn_r(1)/cn_r(5)),x(1)*(mineral_n-thres1))) ;  % above metabolic to active soil
% thres2=x(2)*cn_r(2)/cn_r(5)*ini_litt;
% a_matrix(5,2) = emax * min(1,power(max(1,cn_r(2)/cn_r(5)),x(1)*(mineral_n-thres2))) ;
% thres3=x(2)*cn_r(3)/cn_r(5)*ini_litt;
% a_matrix(5,3) = emax * min(1,power(max(1,cn_r(3)/cn_r(5)),x(1)*(mineral_n-thres3)))*(1- lignin_struc_cmatrix(1)) ;  % above structural to active soil
% thres4=x(2)*cn_r(3)/cn_r(6)*ini_litt;
% a_matrix(6,3) = emax * min(1,power(max(1,cn_r(3)/cn_r(6)),x(1)*(mineral_n-thres4)))*lignin_struc_cmatrix(1); % above structural to slow soil
% thres5=x(2)*cn_r(4)/cn_r(5)*ini_litt;
% a_matrix(5,4) = emax * min(1,power(max(1,cn_r(4)/cn_r(5)),x(1)*(mineral_n-thres5)))*(1 - lignin_struc_cmatrix(2)) ; % below structural to active soil
% thres6=x(2)*cn_r(4)/cn_r(6)*ini_litt;
% a_matrix(6,4) = emax * min(1,power(max(1,cn_r(4)/cn_r(6)),x(1)*(mineral_n-thres6)))*lignin_struc_cmatrix(2); % below structural to slow soil
% %%%%%  Formula 3  %%%%%
thres1=n_cue(1);
a_matrix(5,1) = emax * min(1,power(max(1,cn_r(1)/cn_r(5)),0.54*thres1));  % above metabolic to active soil
thres2=n_cue(2);
a_matrix(5,2) = emax * min(1,power(max(1,cn_r(2)/cn_r(5)),0.54*thres2));
% thres3=x(2)*cn_r(3)/cn_r(5)*ini_litt;
% a_matrix(5,3) = emax * min(1,power(max(1,cn_r(3)/cn_r(5)),x(1)*(mineral_n/thres3-1)))*(1- lignin_struc_cmatrix(1)) ;  % above structural to active soil
% thres4=x(2)*cn_r(3)/cn_r(6)*ini_litt;
% a_matrix(6,3) = emax * min(1,power(max(1,cn_r(3)/cn_r(6)),x(1)*(mineral_n/thres4-1)))*lignin_struc_cmatrix(1); % above structural to slow soil
% thres5=x(2)*cn_r(4)/cn_r(5)*ini_litt;
% a_matrix(5,4) = emax * min(1,power(max(1,cn_r(4)/cn_r(5)),x(1)*(mineral_n/thres5-1)))*(1 - lignin_struc_cmatrix(2)) ; % below structural to active soil
% thres6=x(2)*cn_r(4)/cn_r(6)*ini_litt;
% a_matrix(6,4) = emax * min(1,power(max(1,cn_r(4)/cn_r(6)),x(1)*(mineral_n/thres6-1)))*lignin_struc_cmatrix(2); % below structural to slow soil

% Fixed fraction
% a_matrix(5,1) = frac_soil_metab_aa ;    % above metabolic to active soil
% a_matrix(5,2) = frac_soil_metab_ab ;    % below metabolic to active soil
frac_soil_struct_aa=emax * min(1,power(max(1,cn_r(3)/cn_r(5)),0.54*thres1));
frac_soil_struct_sa=emax * min(1,power(max(1,cn_r(3)/cn_r(6)),0.54*thres1));
a_matrix(5,3) = frac_soil_struct_aa  * (1- lignin_struc_cmatrix(1)) ;  % above structural to active soil
a_matrix(5,4) = frac_soil_struct_ab  * (1 - lignin_struc_cmatrix(2)) ; % below structural to active soil
a_matrix(6,3) = frac_soil_struct_sa * lignin_struc_cmatrix(1); % above structural to slow soil
a_matrix(6,4) = frac_soil_struct_sb * lignin_struc_cmatrix(2); % below structural to slow soil
%
a_matrix(7,5) = frac_passive_active;  % active to passive
a_matrix(6,5) = 1.0 - (0.85-0.68*clay) - a_matrix(7,5);  % active to slow
a_matrix(5,6) = frac_active_slow;  % slow to active
a_matrix(7,6) = frac_passive_slow;  % slow to passive
a_matrix(5,7) = frac_active_passive;  % passive to active
a_matrix(6,7) = frac_slow_passive ;  % passive to slow
% Calibrate CUE when spin-up is over
% N flux between two C/litter pools
nflux_ma(5,1) = a_matrix(5,1)*min(1/cn_r(5),1/cn_r(1)); % above metabolic to active soil
nflux_ma(5,2) = a_matrix(5,2)*min(1/cn_r(5),1/cn_r(2)); % below metabolic to active soil
nflux_ma(5,3) = a_matrix(5,3)*min(1/cn_r(5),1/cn_r(3)); % above structural to active soil
nflux_ma(5,4) = a_matrix(5,4)*min(1/cn_r(5),1/cn_r(4)); % below structural to active soil
nflux_ma(6,3) = a_matrix(6,3)*min(1/cn_r(6),1/cn_r(3)); % above structural to slow soil
nflux_ma(6,4) = a_matrix(6,4)*min(1/cn_r(6),1/cn_r(4)); % below structural to slow soil
nflux_ma(7,5) = a_matrix(7,5)*min(1/cn_r(7),1/cn_r(5)); % active to passive
nflux_ma(6,5) = a_matrix(6,5)*min(1/cn_r(6),1/cn_r(5)); % active to slow
nflux_ma(5,6) = a_matrix(5,6)*min(1/cn_r(5),1/cn_r(6)); % slow to active
nflux_ma(7,6) = a_matrix(7,6)*min(1/cn_r(7),1/cn_r(6)); % slow to passive
nflux_ma(5,7) = a_matrix(5,7)*min(1/cn_r(5),1/cn_r(7)); % passive to active
nflux_ma(6,7) = a_matrix(6,7)*min(1/cn_r(6),1/cn_r(7)); % passive to slow
% Mineralized or immoblized N
for j=1:np
    nmi_ma(j,j)=sum(a_matrix(:,j)/cn_r(j));
end
nmi_ma(5,1)=a_matrix(5,1)*(1/cn_r(5)-1/cn_r(1)); % above metabolic to active soil
nmi_ma(5,2)=a_matrix(5,2)*(1/cn_r(5)-1/cn_r(2)); % above metabolic to active soil
nmi_ma(5,3)=a_matrix(5,3)*(1/cn_r(5)-1/cn_r(3)); % above metabolic to active soil
nmi_ma(5,4)=a_matrix(5,4)*(1/cn_r(5)-1/cn_r(4)); % above metabolic to active soil
nmi_ma(6,3)=a_matrix(6,3)*(1/cn_r(6)-1/cn_r(3)); % above metabolic to active soil
nmi_ma(6,4)=a_matrix(6,4)*(1/cn_r(6)-1/cn_r(4)); % above metabolic to active soil
nmi_ma(7,5)=a_matrix(7,5)*(1/cn_r(7)-1/cn_r(5)); % above metabolic to active soil
nmi_ma(6,5)=a_matrix(6,5)*(1/cn_r(6)-1/cn_r(5)); % above metabolic to active soil
nmi_ma(5,6)=a_matrix(5,6)*(1/cn_r(5)-1/cn_r(6)); % above metabolic to active soil
nmi_ma(7,6)=a_matrix(7,6)*(1/cn_r(7)-1/cn_r(6)); % above metabolic to active soil
nmi_ma(5,7)=a_matrix(5,7)*(1/cn_r(5)-1/cn_r(7)); % above metabolic to active soil
nmi_ma(6,7)=a_matrix(6,7)*(1/cn_r(6)-1/cn_r(7)); % above metabolic to active soil

a_out = a_matrix;

end

     