%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix representation of Century, 7 pools;
% aboveground metabolic litter; belowground meta litter; above structure
% litter; below structure litter; active SOC; slow SOC; passive SOC
% Yuanyuan Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
global np one_year one_day emax
np=7;
one_year=365;
one_day=86400;
emax =0.8;
%========================================================
global clay frac_soil_metab_aa frac_soil_metab_ab frac_soil_struct_aa frac_soil_struct_ab...
    frac_soil_struct_sa frac_soil_struct_sb lignin_struc_cmatrix
global frac_passive_active frac_active_slow frac_passive_slow frac_active_passive frac_slow_passive
global tau_metabolic tau_struct tau_active tau_slow tau_passive flux_tot_coeff
global litter_struct_coef
global cn_r


clay = 0.2;  % 0.2;

frac_soil_metab_aa = 0.45; % aboveground metabolic to active SOC
frac_soil_metab_ab = 0.45; % below metabolic to active SOC
frac_soil_struct_aa =0.55; % above structure to active SOC
frac_soil_struct_ab =0.45; % below structure to active SOC
frac_soil_struct_sa =0.7;  % above structure to slow SOC
frac_soil_struct_sb =0.7;  % below structure to slow SOC

frac_passive_active = 0.004;  % active to passive
frac_active_slow = 0.42;      % slow to active
frac_passive_slow = 0.03;     % slow to passive
frac_active_passive = 0.45;   % passive to active
frac_slow_passive = 0.00;     % passive to slow

lignin_struc_cmatrix(1) = 0.76;  % aboveground lignin in matrix_nextlitter
lignin_struc_cmatrix(2) = 0.72;  % belowground lignin in structure litter

tau_metabolic = 3.5; %0.066*one_year ; %  turnover time in days, calibrated parameters
% tau_metabolic = x(3) ;  % turnover time in days
tau_struct = 0.245*one_year ;
tau_active = 0.149 * one_year;
tau_slow = 5.48* one_year;
tau_passive = 241 * one_year;

flux_tot_coeff = [1.2, 1.4, 0.75];  % only the third is used
litter_struct_coef = 3. ;
cn_r=[30. 30. 100. 100. 12 16 10]; % C/N ratio for [abo_meta,bel_meta,abo_stru,bel_stru,active,slow,passive]

dt=1;
sp=0;  % Spin-up time
ns=200;  %Time length of simulation
soil_weight=1. ;% kg, weight of soil
f=0;
m_fmeta=0.5;
m_fdecomp=296.8; % Parameter for f(Nmin) function

tem=25;
swc=0.6;
matrix_cpools(1,1:7)=[0,0,0,0,0,0,0]; % initial C pool
control_N=[1 1 1 1 1 1 1]; % controlling factor of N on litter decomposition
ini_litt=10; % initial total litter concentration gC/kg soil
ini_soc=10; % initial total soil C concentration gC/kg soil
cn_r(1)=60; % C/N of above-metabolic litter
cn_r(2)=60; % C/N of below-metabolic litter
cn_r(3)=90; % C/N of above-structural litter
cn_r(4)=90; % C/N of above-structural litter
cn_r(5)=12; % C/N of active C pool
cn_r(6)=20; % C/N of slow C pool
cn_r(7)=10; % C/N of passive C pool
mineral_n=0.005;% Concentration of mineral N , g /kg soil
miner_n=mineral_n*soil_weight; % Initial value of mineral N , g
lc_r=0.3; % litter lignin:C ratio
metab_f=max(0.2,0.85-m_fmeta*lc_r);
lignin_struc_cmatrix(1)=min(1,lc_r/(1-metab_f));
lignin_struc_cmatrix(2)=lignin_struc_cmatrix(1);
% fraction of leaf litter belonging to metabolic litter pool
matrix_cpools(1,1)=ini_litt*metab_f;
matrix_cpools(1,3)=ini_litt*(1-metab_f);
matrix_cpools(1,5)=1;
matrix_cpools(1,6)=2;
matrix_cpools(1,7)=7;
matrix_npools(1,1:7)=matrix_cpools(1,1:7)./cn_r(1:7); % initial N pool size
matrix_cpoolsoc(1,1:7)=matrix_cpools(1,1:7);
matrix_cpoolsoc(1,1:4)=[0,0,0,0];
litt_soc_resp(1,1:2)=[0,0];
mineralN(1)=mineral_n; % Soil mineral N concentration (g N/ kg soil)
mat_c_flux=zeros(np,np);   % C flux between different C/litter pools
mat_n_flux=zeros(np,np);      % N flux betweeen different C/litter pools
mat_n_minimb=zeros(np,np);    % mineralization or immobilization of N in litter/carbon pools;

% Start simulation
for ts =1:sp+ns
    matrix_current(1:7) = matrix_cpools(ts,1:7);
    matrix_cursoc(1:7)=matrix_cpoolsoc(ts,1:7);
    
    kk_ma = kk_matrix(ts,sp,tem,swc);  % decomposition x temp scalar x moisture scalar
    n_cue(1)=mineral_n-0.50;
    n_cue(2)=mineral_n-0.50;
    [a_ma,nflux_ma,nmi_ma] = a_matrix(ts,sp,matrix_current,ini_litt,mineral_n,n_cue);  % transfer matrix

    if (a_ma(5,1)/cn_r(5)-1/cn_r(1)<=0)
        control_N(1)=1;
    else
        control_N(1)=min(1,m_fdecomp*mineral_n);
    end
    %
    if (a_ma(5,2)/cn_r(5)-1/cn_r(2)<=0)
        control_N(2)=1;
    else
        control_N(2)=min(1,m_fdecomp*mineral_n);
    end
    %
    n_demand=((1-lignin_struc_cmatrix(1))*a_ma(5,3)/cn_r(5)+lignin_struc_cmatrix(1)*a_ma(6,3)/cn_r(6));
    if (n_demand-1/cn_r(3)<=0)
        control_N(3)=1;
    else
        control_N(3)=min(1,m_fdecomp*mineral_n);
    end
    
    kk_ma(1,1)=kk_ma(1,1)*control_N(1);
    kk_ma(2,2)=kk_ma(2,2)*control_N(2);
    kk_ma(3,3)=kk_ma(3,3)*control_N(3);
    
    matrix_next = matrix_current' + a_ma*kk_ma*matrix_current'*dt; % kk_ma*matrix_current*dt is the decomposition amount;
    matrix_nextsoc=matrix_cursoc' + a_ma*kk_ma*matrix_cursoc'*dt;
    for ipool=1:np
        mat_c_flux(ipool,:)=a_ma(ipool,:).*(kk_ma*matrix_current'*dt)';
        mat_n_flux(ipool,:)=nflux_ma(ipool,:).*(kk_ma*matrix_current'*dt)';
        mat_n_minimb(ipool,:)=nmi_ma(ipool,:).*(kk_ma*matrix_current'*dt)';
    end % for ipool=1:np
    netn=sum(sum(mat_n_minimb)); % negative-mineral, positive-immob
    
    iday=ts*dt;
    tnext = ts +1;
    
    matrix_cpools(tnext,1:7) = matrix_next(1:7,1);
    matrix_cpoolsoc(tnext,1:7)=matrix_nextsoc(1:7,1);
    
    litt_soc_resp(ts,2)=sum(matrix_cpoolsoc(1,1:7))-sum(matrix_cpoolsoc(tnext,1:7));
    litt_soc_resp(ts,1)=sum(matrix_cpools(1,1:7))-sum(matrix_cpools(tnext,1:7))- ...
        litt_soc_resp(ts,2);
    mineralN(tnext)=mineral_n; % Soil mineral N concentration (g N/ kg soil)
    disp(ts);
end %for ts =1:sp+ns

