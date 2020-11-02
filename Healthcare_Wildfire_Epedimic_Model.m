clear; clc;
% simulate the functionality of Healthcare system in Butte county, CA
% After the compounded effect of Wildfire and Epidemic

% Developed by Hassan and Mahmoud 2020 as part of the study
% "Orchestrating Performance of Healthcare Networks Subjected to 
% the Compound Events of Natural Disasters and Pandemic"

% This code is used to run a sample data set for wildfire and pandemic 
% occurred at the time and can't be shared without written permission 
% from authors 

% Model components---------------------------------------------------------
% Hospital 1 = Enloe M.C. % Hospital 2 = Feather river H.
% Hospital 3 = Orchard H. % Hospital 4 = Oroville H.
% Model input--------------------------------------------------------------
% Socio-technical data extracted from: 
% 1- Healthcare providers annual reports prior 2018
% 2- 2018 Census reports
% 3- Camp fire data
% 4- COVID-19 disease transmission data based on SEIR model
% 5- Other soci-technical data from OnTheMap website, Google maps, OSHPD, 
% California Department of Health, insurance, etc.
% Model output-------------------------------------------------------------
% Patient data
% Healthcare system functionality indicators
%%----------------------------------------------------------------------%%

% Model general input 
load('model_input.mat');
% Rating: Hospital rating
% Treat_ability: Ability of the reciever hospital to treat the patients _ assumed based on hospital size 
% Treat_A: Ability of the hospital to treat the patients
% Ambulance: Ambulance of reciever hospital availability
% Air_trans: Air transportation availability 
% Agreement: excisting of hospital agreement
% Pop: Population per census tract
% Private_trans: Private transportation avaialbility
% Insurance: Insurance for each facility
% B_0: Total numer of staffed beds 
% B_em: Emergency rooms 
% B_ICU: Intensive care beds
% B_V: Number of ventilators
% B_in: Inpatient beds 
% B_adm: Admission beds
% staff: Healthcare staff numbers
% TR: Travel time data
% Hospital_travel_time: Travel time between the hospitals
% Housing_DS: Mean social losses for the whole community
% Water: Water recovery
% Power: Power recovery
% Trans: Tranporatation recovery
% Telecom: Telecom recovery
% Westewater: Wastewater recovery
% Fuel: Fuel recovery 

% Total daily patient demand for each census tract and each age group  
load('demand.mat');
% Q_t: number of quarantined cases        
% H_t: number of hospitalized cases
% ICU_t: number of ICU cases
% V_t: number M. ventilators cases
% demand_em: ER demand
% demand_in: Inpatients demand
% demand_ICU: ICU demand
% demand_V: M. ventilators demand

%%----------------------------------------------------------------------%%
%--------------------- Patient-driven model------------------------------%
%%----------------------------------------------------------------------%%
% input and initial values
load('PDM_input.mat','P1','P2','P3_4_5','P6_7','P8','P9','P10','P11','P12');
% P1: Case criticality
% P2: Insurance availability
% P3: Media effect “word of mouth”
% P4: Positive previous experience availability
% P5: Brand availability
% P6: Transportation functionality
% P7: Detour functionality
% P8: Less travel time availability
% P9: Less waiting time
% P10: Hospital ability to treat specific injuries
% P11: Ambulance availability
% P12: Private transportation availability

% Initial value for the components of patient-driven model 
P1_2 = zeros(size(TR,1),length(Ti));  P1_5 = zeros(size(TR,1),length(Ti)); 
P6_8 = zeros(size(TR,1),length(Ti));  P11_12 = zeros(size(TR,1),length(Ti)); 
P9_12 = zeros(size(TR,1),length(Ti)); P_N_t = zeros(size(TR,1),length(Ti)+1);

%%-------------------------------------------------------------------------
% Day 1 patient distribution 
%%-------------------------------------------------------------------------
% Patient probability tree
for i=1:size(P_N_t,1)
P1_2(i,1) = 1- ((1-P1(i,1))*(1-P2(i,1)));      
P1_5(i,1) = P1_2(i,1)*P3_4_5(i,1);
P6_8(i,1) = P6_7(i,1)*P8(i,1);
P11_12(i,1) = 1- ((1-P11(i,1))*(1-P12(i,1)));
P9_12(i,1) = P9(i,1)*P10(i,1)*P11_12(i,1);

P_N_t(i,1) = P1_5(i,1)*P6_8(i,1)*P9_12(i,1);
end

% Patient distribution
load('Patient_dist0.mat','N_em_mod','N_in_mod','N_ICU_mod','N_V_mod',...
    'N_em','N_in','N_ICU','N_V','N_adm','W_t_min','T_t_0','T_t_min');

%%----------------------------------------------------------------------%%
%---------------- Hospital interaction model-----------------------------%
%%----------------------------------------------------------------------%%
% input and initial values
load('HIM_input.mat','I1','I2','I3_4','I5','I6','I7','I8','I9','I10','I11','H_TR');
%I1: Case criticality
%I2: Insurance availability
%I3_4: Transportation
%I5: Travel time
%I6: Agreement
%I7: Telecommunication
%I8: Waiting time
%I9: Treatment ability
%I10: Ambulance
%I11: Air transportation
%H_TR: Hospital travel time

I1_2 = zeros(size(H_TR,1),length(Ti)); 
I3_7 = zeros(size(H_TR,1),length(Ti));
I10_11 = zeros(size(H_TR,1),length(Ti)); 
I8_11 = zeros(size(H_TR,1),length(Ti)); 
I_N_t = zeros(size(H_TR,1),length(Ti));    % Interaction matrix

%%----------------------------------------------------------------------%%
%------- Functionality Fault tree analysis for a hospital cluster--------%
%%----------------------------------------------------------------------%%
 
% ER department quantity functionality input
load('functionality_input_em.mat','R1_em','R2_em','R3_em','R4_em','R5_em','R6_em',...
    'R7_em','R8_em','R9_em','R10_em','R11_em','R12_em','R13_em','R14_em','R15_em',...
    'R16_em','R17_em','R18_em','R19_em','R20_em','R21_em','R22_em','R23_em','R24_em',...
    'R25_em','R26_em','R27_em','R28_em','R29_em','R30_em');
% Inpatient beds quantity functionality input
load('functionality_input_in.mat','R1_in','R2_in','R3_in','R4_in','R5_in','R6_in',...
    'R7_in','R8_in','R9_in','R10_in','R11_in','R12_in','R13_in','R14_in','R15_in',...
    'R16_in','R17_in','R18_in','R19_in','R20_in','R21_in','R22_in','R23_in','R24_in',...
    'R25_in','R26_in','R27_in','R28_in','R29_in','R30_in');
% ICU beds quantity functionality input
load('functionality_input_ICU.mat','R1_ICU','R2_ICU','R3_ICU','R4_ICU','R5_ICU','R6_ICU',...
    'R7_ICU','R8_ICU','R9_ICU','R10_ICU','R11_ICU','R12_ICU','R13_ICU','R14_ICU','R15_ICU',...
    'R16_ICU','R17_ICU','R18_ICU','R19_ICU','R20_ICU','R21_ICU','R22_ICU','R23_ICU','R24_ICU',...
    'R25_ICU','R26_ICU','R27_ICU','R28_ICU','R29_ICU','R30_ICU');
% M. ventilator beds quantity functionality input
load('functionality_input_V.mat','R1_V','R2_V','R3_V','R4_V','R5_V','R6_V',...
    'R7_V','R8_V','R9_V','R10_V','R11_V','R12_V','R13_V','R14_V','R15_V',...
    'R16_V','R17_V','R18_V','R19_V','R20_V','R21_V','R22_V','R23_V','R24_V',...
    'R25_V','R26_V','R27_V','R28_V','R29_V','R30_V');
% R1: Physicians availability
% R2: Nurses availability
% R3: Supporting staff availability
% R4: Alternative staffing availability
% R5: Corridors functionality
% R6: Elevator functionality
% R7: Stairs functionality
% R8: Municipal water functionality
% R9: Backup water system functionality
% R10: Municipal power functionality
% R11: Backup power system functionality
% R12: Transportation network functionality
% R13: Transportation detours availability
% R14: Ambulance service functionality
% R15: Telecommunication service functionality
% R16: Backup telecom service functionality
% R17: Municipal wastewater functionality  
% R18: Backup wastewater functionality
% R19: Drinking water system functionality
% R20: Backup drinking water functionality
% R21: Structural components functionality 
% R22: Non-structural components functionality
% R23: Contents functionality 
% R24: Backup space functionality
% R25: Oxygen availability
% R26: Surgical supply availability
% R27: Rx availability
% R28: Fuel supply availability
% R29: Food supply availability
% R30: Other supply availability


% Patinet-driven model input
load('PDM_mod.mat','P9','P10_em','P10_in','P10_ICU','P10_V');
%P9: Waiting time
%P10_em,P10_in,P10_ICU,P10_V: Different hospital bed availability

% Hospitals interaction model input
load('HIM_mod.mat','I2','I8');
%I2: Insurrance
%I8: Waiting time

for s=1:length(Ti)
for m=1:length(T0) 
 %%----------------------------------------------------------------------%%
% Calculate the ER department quantity functionality
[functionality_c] = functionality_tree(R1_em(m,s),R2_em(m,s),R3_em(m,s),R4_em(m,s),...
    R5_em(m,s),R6_em(m,s),R7_em(m,s),R8_em(m,s),R9_em(m,s),R10_em(m,s),R11_em(m,s),...
    R12_em(m,s),R13_em(m,s),R14_em(m,s),R15_em(m,s),R16_em(m,s),R17_em(m,s),R18_em(m,s),...
    R19_em(m,s),R20_em(m,s),R21_em(m,s),R22_em(m,s),R23_em(m,s),R24_em(m,s),R25_em(m,s),...
    R26_em(m,s),R27_em(m,s),R28_em(m,s),R29_em(m,s),R30_em(m,s));
functionality_em(m,s) =  functionality_c;
B_t_em(m,s) = functionality_em(m,s)*B_em(m);    % total number of available ER staffed beds
 %%----------------------------------------------------------------------%%
% Calculate the inpatient quantity functionality
[functionality_c] = functionality_tree(R1_in(m,s),R2_in(m,s),R3_in(m,s),R4_in(m,s),...
    R5_in(m,s),R6_in(m,s),R7_in(m,s),R8_in(m,s),R9_in(m,s),R10_in(m,s),R11_in(m,s),...
    R12_in(m,s),R13_in(m,s),R14_in(m,s),R15_in(m,s),R16_in(m,s),R17_in(m,s),R18_in(m,s),...
    R19_in(m,s),R20_in(m,s),R21_in(m,s),R22_in(m,s),R23_in(m,s),R24_in(m,s),R25_in(m,s),...
    R26_in(m,s),R27_in(m,s),R28_in(m,s),R29_in(m,s),R30_in(m,s));
functionality_in(m,s) =  functionality_c;
B_t_in(m,s) = functionality_in(m,s)*B_in(m);    % total number of available inpatiet staffed beds
 %%----------------------------------------------------------------------%%
% Calculate the ICU beds quantity functionality
[functionality_c] = functionality_tree(R1_ICU(m,s),R2_ICU(m,s),R3_ICU(m,s),R4_ICU(m,s),...
    R5_ICU(m,s),R6_ICU(m,s),R7_ICU(m,s),R8_ICU(m,s),R9_ICU(m,s),R10_ICU(m,s),R11_ICU(m,s),...
    R12_ICU(m,s),R13_ICU(m,s),R14_ICU(m,s),R15_ICU(m,s),R16_ICU(m,s),R17_ICU(m,s),R18_ICU(m,s),...
    R19_ICU(m,s),R20_ICU(m,s),R21_ICU(m,s),R22_ICU(m,s),R23_ICU(m,s),R24_ICU(m,s),R25_ICU(m,s),...
    R26_ICU(m,s),R27_ICU(m,s),R28_ICU(m,s),R29_ICU(m,s),R30_ICU(m,s));
functionality_ICU(m,s) =  functionality_c;
B_t_ICU(m,s) = functionality_ICU(m,s)*B_ICU(m); % total number of available ICU staffed beds
 %%----------------------------------------------------------------------%%
% Calculate the M. ventilator quantity functionality
[functionality_c] = functionality_tree(R1_V(m,s),R2_V(m,s),R3_V(m,s),R4_V(m,s),...
    R5_V(m,s),R6_V(m,s),R7_V(m,s),R8_V(m,s),R9_V(m,s),R10_V(m,s),R11_V(m,s),...
    R12_V(m,s),R13_V(m,s),R14_V(m,s),R15_V(m,s),R16_V(m,s),R17_V(m,s),R18_V(m,s),...
    R19_V(m,s),R20_V(m,s),R21_V(m,s),R22_V(m,s),R23_V(m,s),R24_V(m,s),R25_V(m,s),...
    R26_V(m,s),R27_V(m,s),R28_V(m,s),R29_V(m,s),R30_V(m,s));
functionality_V(m,s) =  functionality_c;
B_t_V(m,s) = functionality_V(m,s)*B_V(m);       % total number of available M. ventilator staffed beds
%%----------------------------------------------------------------------%%
 % Total quantitative functionality %Hospital functionality
 functionality(m,s) = (functionality_in(m,s)*B_in(m)+functionality_em(m,s)*B_em(m)...
                     +functionality_ICU(m,s)*B_ICU(m)+functionality_V(m,s)*B_V(m))...
                     /(B_in(m)+B_em(m)+B_ICU(m)+B_V(m));
%%----------------------------------------------------------------------%%               
if s==1
     kk = 1.0;
 else
 kk = find(hospital_census(:,s-1)==(m)); % to get the census tracts anticipated to each hospital
 end
 H_Travel_time(m,s) = mean(TR(kk,s));    % to get the travel time associated with each hospital
 H_Travel_time(isnan(H_Travel_time)) = mean_H_Travel_time; 
 clear kk
%%-----------------------------------------------------------------------%%
% Qualitative functionality
% Waiting time
W_t(m,s) = max(a_0(m)+H_Travel_time(m,s)+a_t*(B_em(m)-B_t_em(m,s))/B_em(m)+a_e*(N_em(m,s)-N_0_em(m))/N_em(m),a_0(m)); %waiting time in ER departments
% Treatment time
T_t(m,s) = (R1_em(m,s)/N_em(m,s));
% Hospitalization service accessibility
S_A(m,s) = min(max((W_max - W_t(m,s))/(W_max - W_t_min(m)),0.0),1.0);
% Hospitalization service effectiveness
S_E(m,s) = min(max((T_t(m,s) - T_t_min)/(T_t_0(m) - T_t_min),0.0),1.0);
% Qualitative functionality 
Q_s(m,s) = S_A(m,s)*S_E(m,s);
% Total functionality 
F(m,s) = Q_s(m,s).^(I_quality(s))*functionality(m,s);      
end

%%----------------------------------------------------------------------%%
%--------------------- Patient-driven model------------------------------%
%%----------------------------------------------------------------------%%

for i=1:size(P_N_t,1)
P1_2(i,s) = 1- ((1-P1(i,s))*(1-P2(i,s)));      
P1_5(i,s) = P1_2(i,s)*P3_4_5(i,s);
P6_8(i,s) = P6_7(i,s)*P8(i,s);
P11_12(i,s) = 1- ((1-P11(i,s))*(1-P12(i,s)));
P9_12(i,s) = P9(i,s)*P11_12(i,s);

P_N_em(i,s) = P1_5(i,s)*P6_8(i,s)*P9_12(i,s)*P10_em(i,s);
P_N_in(i,s) = P1_5(i,s)*P6_8(i,s)*P9_12(i,s)*P10_in(i,s);
P_N_ICU(i,s) = P1_5(i,s)*P6_8(i,s)*P9_12(i,s)*P10_ICU(i,s);
P_N_V(i,s) = P1_5(i,s)*P6_8(i,s)*P9_12(i,s)*P10_V(i,s);
end

[N_em_mod_c,N_in_mod_c,N_ICU_mod_c,N_V_mod_c,N_em_c,N_in_c,N_ICU_c,N_V_c,N_adm_c,hospital_census_c]=...
    Patient_dist(T0,s,Pop,P_N_em,P_N_in,P_N_ICU,P_N_V,demand_em,demand_in,demand_ICU,demand_V);
    hospital_census(:,s) = hospital_census_c;
    N_em(:,s+1)  = N_em_c;
    N_in(:,s+1)  = N_in_c;
    N_ICU(:,s+1) = N_ICU_c;
    N_V(:,s+1)   = N_V_c;
    N_adm(:,s+1) = N_adm_c;

%%----------------------------------------------------------------------%%
%---------------- Hospital interaction model-----------------------------%
%%----------------------------------------------------------------------%%

 for i=1:size(H_TR,1)*2
 I1_2(i,s) = 1- ((1-I1(i,s))*(1-I2(i,s)));
 I3_7(i,s) = I3_4(i,s)*I5(i,s)*I6(i,s)*I7(i,s);
 I10_11(i,s) = 1- ((1-I10(i,s))*(1-I11(i,s)));
 I8_11(i,s) = I10_11(i,s)*I8(i,s)*I9(i,s);
 
 I_N_t(i,s) = I1_2(i,s)*I3_7(i,s)*I8_11(i,s);

 end

Interaction = [0           I_N_t(1,s)  I_N_t(2,s)  I_N_t(3,s);...
               I_N_t(7,s)  0           I_N_t(4,s)  I_N_t(5,s);...
               I_N_t(8,s)  I_N_t(9,s)  0           I_N_t(6,s);...
               I_N_t(10,s) I_N_t(11,s) I_N_t(12,s) 0]; % Interaction matrix

% Patient redistribution to consider the hospital interaction
[N_in_c,N_em_c,N_ICU_c,N_V_c,N_adm_c] = HIM_redist(s,Interaction,N_in,tr_in,N_0_in,B_t_in,...
    N_em,tr_em,N_0_em,B_t_em,N_ICU,tr_ICU,N_0_ICU,B_t_ICU,N_V,tr_V,N_0_V,B_t_V);
N_in(:,s+1) = N_in_c(:,s+1);
N_em(:,s+1) = N_em_c(:,s+1);
N_ICU(:,s+1) = N_ICU_c(:,s+1);
N_V(:,s+1) = N_V_c(:,s+1);
N_adm(:,s+1) = N_adm_c(:,s+1);

% Patient overflow calculations
[N_in_overflow_c,N_em_overflow_c,N_ICU_overflow_c,N_V_overflow_c,N_adm_overflow_c]...
    = Overflow(s,T0,N_in,tr_in,B_t_in,N_em,tr_em,Exp_em,B_t_em,N_ICU,tr_ICU,B_t_ICU,N_V,tr_V,B_t_V);
N_in_overflow (:,s+1)  = N_in_overflow_c(:,s+1); 
N_em_overflow (:,s+1)  = N_em_overflow_c(:,s+1);
N_ICU_overflow (:,s+1) = N_ICU_overflow_c(:,s+1);
N_V_overflow (:,s+1)   = N_V_overflow_c(:,s+1);
N_adm_overflow (:,s+1) = N_adm_overflow_c(:,s+1);

% Total patients overflow
Overflow_in(s+1)  = sum(N_in_overflow(:,s+1));
Overflow_em(s+1)  = sum(N_em_overflow(:,s+1));
Overflow_ICU(s+1) = sum(N_ICU_overflow(:,s+1));
Overflow_V(s+1)   = sum(N_V_overflow(:,s+1));
Overflow_adm(s+1) = sum(N_adm_overflow(:,s+1));

% To adjust the assummed patient number at the first day 
if s==1
 N_in(:,1) = N_in(:,2); N_em(:,1) = N_em(:,2); N_ICU(:,1) = N_ICU(:,2); 
 N_V(:,1) = N_V(:,2);  N_adm(:,1) = N_adm(:,2);
end   
end

clear W_t_min
W_t_min = zeros(1,length(T0));
H_Travel_time(:,1) = H_Travel_time(:,2);
for i=1:length(T0)
W_t_min(i) = a_0(i)+H_Travel_time(i,1); % min waiting time
T_t_0(i) = (R1_em(i,1)/N_0_em(i));      % treatment time for each patient (day)
end

for s=1:length(Ti)
for m=1:length(T0)
N_in(m,s)  = N_in(m,s)  - N_in_overflow (m,s);
N_em(m,s)  = N_em(m,s)  - N_em_overflow (m,s);
N_ICU(m,s) = N_ICU(m,s) - N_ICU_overflow (m,s);
N_V(m,s)   = N_V(m,s)   - N_V_overflow (m,s);
N_adm(m,s) = N_adm(m,s) - N_adm_overflow (m,s);
end 
end

for s=1:length(Ti)
for m=1:length(T0)
%%-----------------------------------------------------------------------%%
% Qualitative functionality re-calculations
B_t_em(m,s) = functionality_em(m,s)*B_em(m); % number of available staffed beds
W_t(m,s) = max(a_0(m)+H_Travel_time(m,s)+a_t*(B_em(m)-B_t_em(m,s))/B_em(m)+...
    a_e*(N_em(m,s)+N_em_overflow(m,s)-N_em(m,1))/N_em(m,1),W_t_min(m)); %waiting time
T_t(m,s) = (R1_em(m,s)/N_em(m,s)); % treatment time
S_A(m,s) = min(max((W_max - W_t(m,s))/(W_max - W_t_min(m)),0.0),1.0); % Accessibility
S_E(m,s) = min(max((T_t(m,s) - T_t_min)/(T_t_0(m) - T_t_min),0.0),1.0); % Effectiveness
% Qualitative functionality re-calculations
Q_s(m,s) = S_A(m,s)*S_E(m,s);  
% Total functionality re-calculations
F(m,s) = Q_s(m,s).^(I_quality(s))*functionality(m,s);      
end
end
clear hospital_census functionality_c
% Results collection section
F_H(:) = (F(1,:).*B_0(1)+F(2,:).*B_0(2)+F(3,:).*B_0(3)+F(4,:).*B_0(4))/sum(B_0);         % Hopsital total func.
F_Q(:) = (Q_s(1,:).*B_0(1)+Q_s(2,:).*B_0(2)+Q_s(3,:).*B_0(3)+Q_s(4,:).*B_0(4))/sum(B_0); % Hopsital quality func.
F_W(:) = (W_t(1,:).*B_t_em(1,:)+W_t(2,:).*B_t_em(2,:)+W_t(3,:).*B_t_em(3,:)+W_t(4,:).*B_t_em(4,:))./sum(B_t_em); % Hopsital waiting time
F_T(:) = (T_t(1,:).*B_t_em(1,:)+T_t(2,:).*B_t_em(2,:)+T_t(3,:).*B_t_em(3,:)+T_t(4,:).*B_t_em(4,:))./sum(B_t_em); % Hopsital treatment time

F_N_in(:) = N_in(1,:)+N_in(2,:)+N_in(3,:)+N_in(4,:);       % Total number of inpatients 
F_N_em(:) = N_em(1,:)+N_em(2,:)+N_em(3,:)+N_em(4,:);       % Total number of ER patients
F_N_ICU(:) = N_ICU(1,:)+N_ICU(2,:)+N_ICU(3,:)+N_ICU(4,:);  % Total number of ICU admissions
F_N_V(:) = N_V(1,:)+N_V(2,:)+N_V(3,:)+N_V(4,:);            % Total number of mechanical ventilator patients 
F_N_in_overflow(:) = N_in_overflow(1,:)+N_in_overflow(2,:)+N_in_overflow(3,:)+N_in_overflow(4,:);       % number of untreated inpatients 
F_N_em_overflow(:) = N_em_overflow(1,:)+N_em_overflow(2,:)+N_em_overflow(3,:)+N_em_overflow(4,:);       % number of untreated ER patients 
F_N_ICU_overflow(:) = N_ICU_overflow(1,:)+N_ICU_overflow(2,:)+N_ICU_overflow(3,:)+N_ICU_overflow(4,:);  % number of untreated ICU patients
F_N_V_overflow(:) = N_V_overflow(1,:)+N_V_overflow(2,:)+N_V_overflow(3,:)+N_V_overflow(4,:);            % number of untreated mechanical ventilator patients  

%%-----------------------------------------------------------------------%%
%-------------------------- Plot Section --------------------------------%%
%%-----------------------------------------------------------------------%%
Plot_var = N_in;
Plot_dot = N_in_overflow;
Plot_var_t = sum(Plot_var,1);
Plot_dot_t = sum(Plot_dot,1);

model_1 = plot(0:length(Ti)-1,Plot_var(1,1:length(Ti))+Plot_dot(1,1:length(Ti)),'--r');
hold on
model_2 = plot(0:length(Ti)-1,Plot_var(2,1:length(Ti))+Plot_dot(2,1:length(Ti)),'--g');
model_3 = plot(0:length(Ti)-1,Plot_var(3,1:length(Ti))+Plot_dot(3,1:length(Ti)),'--b');
model_4 = plot(0:length(Ti)-1,Plot_var(4,1:length(Ti))+Plot_dot(4,1:length(Ti)),'--k');
model_5 = plot(0:length(Ti)-1,Plot_var(1,1:length(Ti)),'-r');
model_6 = plot(0:length(Ti)-1,Plot_var(2,1:length(Ti)),'-g');
model_7 =plot(0:length(Ti)-1,Plot_var(3,1:length(Ti)),'-b');
model_8 = plot(0:length(Ti)-1,Plot_var(4,1:length(Ti)),'-k');

set(model_1,'LineStyle','--','LineWidth',1.0,'Color',[1.0 0.0 0.6]);
set(model_5,'LineStyle','-','LineWidth',2.0,'Color',[1.0 0.0 0.6]);
set(model_2,'LineStyle','--','LineWidth',1.0,'Color',[1.0 0.6 0.6]);
set(model_6,'LineStyle','-','LineWidth',2.0,'Color',[1.0 0.6 0.6]);
set(model_3,'LineStyle','--','LineWidth',1.0,'Color',[1.0 0.8 0.6]);
set(model_7,'LineStyle','-','LineWidth',2.0,'Color',[1.0 0.8 0.6]);
set(model_4,'LineStyle','--','LineWidth',1.0,'Color',[1.0 0.4 0.6]);
set(model_8,'LineStyle','-','LineWidth',2.0,'Color',[1.0 0.4 0.6]);
hXLabel = xlabel('Time (day)');
hYLabel = ylabel('Inpatient number');

set( gca                       , ...
    'FontName'   , 'AvantGarde' );
set([ hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde' ,...
    'FontWeight' , 'bold'      );
set([gca]             , ...
    'FontSize'   , 16           );
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 18          );
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:100:600 ,...
  'XTick'       , 0:40:360 ,...
  'FontSize'    , 16     ,...
  'LineWidth'   , 1         );
xlim([0 360])
ylim([0 600])


% Success tree calculations for the hospital functionality 
function [functionality] = functionality_tree(R1,R2,R3,R4,R5,R6,R7,...
    R8,R9,R10,R11,R12,R13,R14,R15,R16,R17,R18,...
    R19,R20,R21,R22,R23,R24,R25,R26,R27,R28,R29,R30)
% Fault tree analysis 
 R1_3 = R1*R2*R3;
 R1_4 = 1- ((1-R1_3)*(1-R4));                    % Staff availability
 
 R6_7 = 1- ((1-R6)*(1-R7));   
 R5_7 = R5*R6_7;                                  % Accecability 
 
 R8_9 = 1- ((1-R8)*(1-R9));
 R10_11 = 1- ((1-R10)*(1-R11));
 R12_13 = 1- ((1-R12)*(1-R13));
 R12_14 = R12_13*R14;
 R15_16 = 1- ((1-R15)*(1-R16));
 R17_18 = 1- ((1-R17)*(1-R18));
 R19_20 = 1- ((1-R19)*(1-R20));
 R8_20 = R8_9*R10_11*R12_14*R15_16*R17_18*R19_20; %Supportive infrastructure
 R21_23 = R21*R22*R23; 
 R21_24 = 1- ((1-R21_23)*(1-R24));                % Working space
 R5_24 = R5_7*R8_20*R21_24;                       % Spcae availability   
 
 R25_30 = R25*R26*R27*R28*R29*R30;                % Supplies availability
 functionality =  R1_4*R5_24*R25_30;              % Quantity functionality
end

% Patinet distribution function based on Patient-driven model data
function [N_em_mod,N_in_mod,N_ICU_mod,N_V_mod,N_em,N_in,N_ICU,N_V,N_adm,hospital_census]=...
    Patient_dist(T0,s,Pop,P_N_em,P_N_in,P_N_ICU,P_N_V,demand_em,demand_in,demand_ICU,demand_V)
for i=1:length(Pop)  
pufer = find(P_N_em(4*i-3:4*i,s)==max(P_N_em(4*i-3:4*i,s))); hospital_census(i) = pufer(1);
end

for i=1:length(T0)
     for j=1:length(Pop)
    N_em_mod(j,i) = (P_N_em(4*j-4+i,s)/sum(P_N_em(4*j-3:4*j,s))*demand_em(j,s+1));
    N_in_mod(j,i) = (P_N_in(4*j-4+i,s)/sum(P_N_in(4*j-3:4*j,s))*demand_in(j,s+1));
    N_ICU_mod(j,i) = (P_N_ICU(4*j-4+i,s)/sum(P_N_ICU(4*j-3:4*j,s))*demand_ICU(j,s+1));
    N_V_mod(j,i) = (P_N_V(4*j-4+i,s)/sum(P_N_V(4*j-3:4*j,s))*demand_V(j,s+1));
     end
    N_em(i) = sum(N_em_mod(:,i));
    N_in(i) = sum(N_in_mod(:,i));
    N_ICU(i) = sum(N_ICU_mod(:,i));
    N_V(i) = sum(N_V_mod(:,i));
    N_adm(i) = N_in(i)+N_ICU(i)+N_V(i);
end
end

% Function to calculate the patient overflow
function [N_in_overflow,N_em_overflow,N_ICU_overflow,N_V_overflow,N_adm_overflow] = Overflow(s,T0,N_in,tr_in,B_t_in,...
    N_em,tr_em,Exp_em,B_t_em,N_ICU,tr_ICU,B_t_ICU,N_V,tr_V,B_t_V)
for i=1:length(T0)
N_in_overflow (i,s+1) = max(sum(N_in(i,s+1))- sum(tr_in(s+1)*B_t_in(i,s)),0); 
N_em_overflow (i,s+1) = max(sum(N_em(i,s+1))- sum(tr_em(s+1)*Exp_em*B_t_em(i,s)),0); 
N_ICU_overflow (i,s+1) = max(sum(N_ICU(i,s+1))- sum(tr_ICU(s+1)*B_t_ICU(i,s)),0); 
N_V_overflow (i,s+1) = max(sum(N_V(i,s+1))- sum(tr_V(s+1)*B_t_V(i,s)),0);   
N_adm_overflow (i,s+1) = N_in_overflow (i,s+1)+N_ICU_overflow (i,s+1)+N_V_overflow (i,s+1);
end
end

% Function to redistribute the patients
function [N_in,N_em,N_ICU,N_V,N_adm] = HIM_redist(s,Interaction,N_in,tr_in,N_0_in,B_t_in,...
    N_em,tr_em,N_0_em,B_t_em,N_ICU,tr_ICU,N_0_ICU,B_t_ICU,N_V,tr_V,N_0_V,B_t_V)

% Inpatient transfer
    Over_capacity_in = find(N_in(:,s+1)>tr_in(s+1)*N_0_in(:)); less_capacity_in = find(N_in(:,s+1)<=1.0*N_0_in(:));
    if isempty(Over_capacity_in) == 0
    for j=Over_capacity_in'
    for i=less_capacity_in'
        less_capacity_puffer = max(min(N_0_in(i)-N_in(i,s+1),tr_in(s+1)*B_t_in(i,s)),0);
        tranfered_in(j,s+1) = max(min(Interaction(j,i)'.*(N_in(j,s+1)-tr_in(s+1)*N_0_in(j))/sum(Interaction(j,less_capacity_in)),less_capacity_puffer),0);
        N_in(i,s+1) = round(min(N_in(i,s+1)+tranfered_in(j,s+1),N_0_in(i)));
        N_in(j,s+1) = round(N_in(j,s+1)-tranfered_in(j,s+1));
    end
    end
    else 
       tranfered_in(:,s+1) = 0.0;
    end
% Emeregency transfer
    Over_capacity_em = find(N_em(:,s+1)>tr_em(s+1)*N_0_em(:)); less_capacity_em = find(N_em(:,s+1)<=1.0*N_0_em(:));
    if isempty(Over_capacity_em) == 0
    for j=Over_capacity_em'
    for i=less_capacity_em'
        less_capacity_puffer = max(min(N_0_em(i)-N_em(i,s+1),tr_em(s+1)*B_t_em(i,s)),0);
        tranfered_em(j,s+1) = max(min(Interaction(j,i)'.*(N_em(j,s+1)-tr_em(s+1)*N_0_em(j))/sum(Interaction(j,less_capacity_em)),less_capacity_puffer),0);
        N_em(i,s+1) = round(min(N_em(i,s+1)+tranfered_em(j,s+1),N_0_em(i)));
        N_em(j,s+1) = round(N_em(j,s+1)-tranfered_em(j,s+1));
    end
    end
    else 
       tranfered_em(:,s+1) = 0.0;
    end
% ICU transfer
    Over_capacity_ICU = find(N_ICU(:,s+1)>tr_ICU(s+1)*N_0_ICU(:)); less_capacity_ICU = find(N_ICU(:,s+1)<=1.0*N_0_ICU(:));
    if isempty(Over_capacity_ICU) == 0
    for j=Over_capacity_ICU'
    for i=less_capacity_ICU'
        less_capacity_puffer = max(min(N_0_ICU(i)-N_ICU(i,s+1),tr_ICU(s+1)*B_t_ICU(i,s)),0);
        tranfered_ICU(j,s+1) = max(min(Interaction(j,i)'.*(N_ICU(j,s+1)-tr_ICU(s+1)*N_0_ICU(j))/sum(Interaction(j,less_capacity_ICU)),less_capacity_puffer),0);
        N_ICU(i,s+1) = round(min(N_ICU(i,s+1)+tranfered_ICU(j,s+1),N_0_ICU(i)));
        N_ICU(j,s+1) = round(N_ICU(j,s+1)-tranfered_ICU(j,s+1));
    end
    end
    else 
       tranfered_ICU(:,s+1) = 0.0;
    end
% Ventilator transfer
    Over_capacity_V = find(N_V(:,s+1)>tr_V(s+1)*N_0_V(:)); less_capacity_V = find(N_V(:,s+1)<=1.0*N_0_V(:));
    if isempty(Over_capacity_V) == 0
    for j=Over_capacity_V'
    for i=less_capacity_V'
        less_capacity_puffer = max(min(N_0_V(i)-N_ICU(i,s+1),tr_V(s+1)*B_t_V(i,s)),0);
        tranfered_V(j,s+1) = max(min(Interaction(j,i)'.*(N_V(j,s+1)-tr_V(s+1)*N_0_V(j))/sum(Interaction(j,less_capacity_V)),less_capacity_puffer),0);
        N_V(i,s+1) = round(min(N_V(i,s+1)+tranfered_V(j,s+1),N_0_V(i)));
        N_V(j,s+1) = round(N_V(j,s+1)-tranfered_V(j,s+1));
    end
    end
    else 
       tranfered_V(:,s+1) = 0.0;
    end
    
N_adm(:,s+1) = N_in(:,s+1)+N_ICU(:,s+1)+N_V(:,s+1);

end