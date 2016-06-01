clc
close all
clear all
format long g
%-----------------------Definitions:---------------------------------------

%Rocket information. All units in: [kg, kJ, kPa, Kelvin, degrees, mol]
%Defining: Chamber volume, ambient pressure, ambient temperature, 
%oxidizer purity percentage, oxidizer used, oxidizer mass flow, oxidizer mass, 
%gas constant, chamber area, throat area, calculation timesteps, counter
%enthalpy for: water, O2, CO2, Stefan-Boltzman constant, Avogadros number
%injection duration, gamma
V_chamber   = 3;            %[L]
P_amb       = 101.3;            %[kPa]
T_amb       = 273.15 + 20;      %[K]
x_oxidizer  = 0.8;              %fraction
f           = 0.59;             %fraction
m_dot_oxidizer = 0.246;         %[kg/s]
mass_H2O2   = 1.3 * x_oxidizer; %[kg]
mass_H2O    = 1.3 - mass_H2O2;  %[kg]
R           = 8.3145;           %[kJ/(mol*K)]
A_chamber   = (0.094)^2;        %[m]^2
A_throat    = (0.0213)^2;       %[m]^2
A_exit      = (0.0337)^2;       %[m]^2
dt          = 0.001;             %[s]
k           = 2;                %[s]
H_water     = 83.93;            %[kJ/kg]
H_oxygen    = -4.887;           %[kJ/kg]
H_CO2       = -5.168;           %[kJ/kg]
kB          = 1.38065E-23;      %[J/K]
NA          = 6.02214086E23;    %[mol]^-1
time_state1 = (mass_H2O2+mass_H2O)/m_dot_oxidizer; %[s]
gamma       = 1.2;

%Defining: H2O2, H2O, O2, Plastic grain, CO2 Molar masses:
M_H2O2 = 0.0340147; %[kg/mol]
M_H2O  = 0.0180153; %[kg/mol]
M_O2   = 0.0319988; %[kg/mol]
M_PLA  = 0.0720000; %[kg/mol]
M_CO2  = 0.0440095; %[kg/mol]
M_air  = 0.0289700; %[kg/mol]

%Defining specific enthalpy from decomposition and combustion reaction:
%2*H2O2 -> 2*H2O + O2 and C3H4O2 + 3*O2 -> 3*CO2 + 2*H2O respectively
H2O2_Gibbs_free      = 98.2;                   %[kJ/mol]
DELTAh_decomposition = H2O2_Gibbs_free/M_H2O2; %[kJ/kg]
DELTAh_combustion    = 18000;                  %[kJ/kg]

% ---------------------- Flow rates ---------------------------------------

%-----------------------State 1: Injection:--------------------------------
%Injection starts at t=0
n_H2O2_1    = mass_H2O2/M_H2O2;                    %[mol]
n_H2O_1     = mass_H2O/M_H2O;                      %[mol]
n_O2_1      = 0;
n_CO2_1     = 0;

n_dot_H2O2_1 = x_oxidizer*m_dot_oxidizer/M_H2O2;
n_dot_H2O_1  = (1-x_oxidizer)*m_dot_oxidizer/M_H2O;
n_dot_O2_1   = 0;
n_dot_CO2_1  = 0;


n_flow_1    = n_H2O2_1 + n_H2O_1;                         %Total flow from injection
M_1         = (n_H2O_1*M_H2O + n_H2O2_1*M_H2O2)/n_flow_1; %mass of injected matter
m_1         = (n_H2O_1*M_H2O + n_H2O2_1*M_H2O2);

%-----------------------State 2: Decomposition:----------------------------
%Decomposition starts one time step dt later:

n_H2O2_2 = 0;                   %[mol]
n_H2O_2  = n_H2O_1 +  n_H2O2_1; %[mol]
n_O2_2   = 0.5 * n_H2O2_1;      %[mol]
n_CO2_2  = 0;                   %[mol]

n_dot_H2O2_2 = 0;
n_dot_H2O_2  = n_dot_H2O_1 + n_dot_H2O2_1;
n_dot_O2_2   = 0.5*n_dot_H2O2_1;
n_dot_CO2_2  = 0;

n_flow_2 = n_H2O2_2 + n_H2O_2 + n_O2_2;     %[mol]
m_2      = n_H2O_2 * M_H2O + n_O2_2 * M_O2; %[kg]

%-----------------------State 3: Combustion:-------------------------------
% starts one time step dt later yet:

n_H2O2_3 = 0;
n_H2O_3  = n_H2O_2 +  f*2/3*n_O2_2; %[mol]
n_O2_3   = (1-f)*n_O2_2;
n_CO2_3  = f*n_O2_2;

n_dot_H2O2_3 = 0;
n_dot_H2O_3  = n_dot_H2O_2 + f * 2/3 * n_dot_O2_2;
n_dot_O2_3   = (1-f)*n_dot_O2_2;
n_dot_CO2_3  = f*n_dot_O2_2;

n_flow_3 = n_H2O2_3 + n_H2O_3 + n_O2_3 +n_CO2_3;
m_PLA_3  = n_CO2_3/3*M_PLA;         %Mass of burned material. 
m_dot_PLA_3  = n_dot_CO2_3/3*M_PLA; %Only CO2 and H2O comes from burning
m_3      = n_H2O_3 * M_H2O + n_O2_3 * M_O2 + n_CO2_3 * M_CO2;
M_3      = m_3/n_flow_3;

%----------------------------Before rocket start:--------------------------

P_tot(1) = P_amb;
T_tot(1) = T_amb;
n_tot(1) = (P_tot(1)*V_chamber)/(R*T_tot(1));


%----------------------------Timestep dt later:----------------------------
%Initial injection occurs

n_tot(2) = n_tot(1) + n_flow_1 * dt
T_tot(2) = T_amb;
P_tot(2) = n_tot(2)*R*T_tot(2)/V_chamber;

%----------------------------Timestep dt later:----------------------------
%Initial decomposition occurs

%Enthalpy released in decomposition:
H_ref_2 = n_dot_H2O_2 * M_H2O * H_water * dt + n_dot_O2_2 * M_O2 * H_oxygen * dt + n_dot_CO2_2 * M_CO2 * H_CO2 * dt;
H_2     = H_ref_2 + DELTAh_decomposition * n_dot_H2O2_1 * M_H2O2 * dt;


%---------------------------INITIATE EES TRANSFER--------------------------
tempchan = ddeinit('EES','DDE');
%Opens EES work-file
EESload  = ddeexec(tempchan,'[Open enthalpytotemp2.EES]');

    P_temp = P_tot(2);
    H_temp = H_2;
    n_tempH2O = n_dot_H2O_2 * M_H2O * dt;
    n_tempO2  = n_dot_O2_2 * M_O2 * dt;
    n_tempCO2 = n_dot_CO2_2 * M_CO2 * dt;
    n_tempAir = n_tot(1);

save enthalpytotemp.txt P_temp H_temp n_tempH2O n_tempO2 n_tempCO2 n_tempAir -ascii;
EESload = ddeexec(tempchan,'[Solve]');

    T(2) = csvread('tempin.csv');



%--------------------- Figure Plotting ------------------------------------
total_time = length(P_tot)*dt;
tspan = linspace(0,total_time,length(P_tot));

figure(1)
    yyaxis left
    title('Pressure and temperature over time')
    plot(tspan,P_tot/100,'-o')
    hold on
    yyaxis right 
    plot(tspan,T_tot-T_amb)
    
    axis([0 max(tspan) 0 max(P_tot)/100])
    ylabel('Pressure [bar]')
    xlabel('time [s]')



