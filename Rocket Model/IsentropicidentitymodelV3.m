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
T_auto_MDF  = 273.15 + 220;     %[K]

%Defining: H2O2, H2O, O2, Plastic grain, CO2 Molar masses:
M_H2O2 = 0.0340147; %[kg/mol]
M_H2O  = 0.0180153; %[kg/mol]
M_O2   = 0.0319988; %[kg/mol]
M_PLA  = 0.0720000; %[kg/mol]
M_CO2  = 0.0440095; %[kg/mol]
M_air  = 0.0289700; %[kg/mol]
M_Nit   = 0.0280134; %[kg/mol]

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
m_inside(1) = n_tot(1)*M_air;

m_out(1) = 0;

delta_m(1)  = m_inside(1) - m_out(1)*dt;

%----------------------------Timestep dt later:----------------------------
%Initial injection occurs

n_tot(2) = n_tot(1) + n_flow_1/time_state1 * dt;
T_tot(2) = T_amb;
P_tot(2) = n_tot(2)*R*T_tot(2)/V_chamber;
m_inside(2) = n_tot(1)*M_air+m_1/time_state1 * dt;

m_out(2) = 0;

delta_m(2)  = m_inside(2) - m_out(2)*dt;

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
    n_tempH2O = n_dot_H2O_2 * dt;
    n_tempO2  = n_dot_O2_2 * dt + 0.79 * n_tot(1);
    n_tempCO2 = n_dot_CO2_2 * dt;
    n_tempNit = 0.21 * n_tot(1);

save enthalpytotemp.txt P_temp H_temp n_tempH2O n_tempO2 n_tempCO2 n_tempNit -ascii;
EESload = ddeexec(tempchan,'[Solve]');

    T(3) = csvread('tempin.csv');
ddeterm(tempchan);

%Calculating new amount of substance
n_tot(3) = n_tempH2O + n_tempO2 + n_tempCO2 + n_tempNit;

%calculating working pressure for next EES transfer
P(3) = n_tot(3) * R * T(3)/V_chamber;
H(3) = H_temp;

%Calculating how much matter is blown out of the rocket
T_e(3)   = (P_amb./P(3))^(1-1/gamma) * T(3);
rho_e(3) = P_amb./(R*T_e(3));
v_e(3)   = sqrt(2*(P(3)-P_amb)./rho_e(3));
m_out(3) = A_exit*v_e(3)*rho_e(3);

%calculating mass of matter inside rocket
m_inside(3) = n_tempH2O * M_H2O + n_tempO2 * M_O2 + n_tempCO2 * M_CO2 + n_tempNit * M_Nit;

%finding change in mass from matter flowing out. m_out has units [kg/s],
%thus multiplying by dt is required
delta_m(3)    = m_inside(3) - m_out(3)*dt;
  
%Air leaves the rocket first, as it is closest to the nozzle
n_air = n_tot(1) - m_out(3)/M_air * dt;

%---------------------------INITIATE EES TRANSFER--------------------------
tempchan = ddeinit('EES','DDE');
%Opens EES work-file
EESload  = ddeexec(tempchan,'[Open enthalpytotemp2.EES]');

%Loads in temporary values to save to EES with air abundances reducing per.
%timestep.
    P_temp = P(3);
    n_tempH2O = n_dot_H2O_2 * dt;
    n_tempO2  = n_dot_O2_2 * dt + 0.79 * n_air
    n_tempCO2 = n_dot_CO2_2 * dt;
    n_tempNit = 0.21 * n_air;

    H_tempref = n_tempH2O * M_H2O * H_water + n_tempO2 * M_H2O * H_oxygen + n_tempCO2 * H_CO2; %+ n_tempNit * H_water
    H_temp    = H_tempref + DELTAh_decomposition * n_dot_H2O2_1 * M_H2O2 * dt;
%Saves variables in file enthalpytotemp.txt in ascii format. 
save enthalpytotemp.txt P_temp H_temp n_tempH2O n_tempO2 n_tempCO2 n_tempNit -ascii;

%Loads data into EES and performs the function [solve]
EESload = ddeexec(tempchan,'[Solve]');
%EES loads the data, reads vars, calculates T based on enthalpy, saves
%temperature in tempin.csv which is loaded into T_tot var below.
    T_tot(3) = csvread('tempin.csv');
%Terminates data exchange
ddeterm(tempchan);

H_tot(3) = H_temp
P_tot(3) = n_tot(3)*R*T_tot(3)/V_chamber;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------Initiate for-loop until T=T_auto:-------------------------
k = 3;
while max(T_tot)<T_auto_MDF
    k = k + 1
    %---------------------------INITIATE EES TRANSFER--------------------------
    tempchan = ddeinit('EES','DDE');
%Opens EES work-file
    EESload  = ddeexec(tempchan,'[Open enthalpytotemp2.EES]');

        P_temp = P_tot(k-1);
        H_temp = H_tot(k-1);
        n_tempH2O = n_tempH2O + n_dot_H2O_2 * dt;
        n_tempO2  = n_tempO2 + n_dot_O2_2 * dt;
        n_tempCO2 = n_tempCO2 + n_dot_CO2_2 * dt;
        n_tempNit = n_tempNit;

    save enthalpytotemp.txt P_temp H_temp n_tempH2O n_tempO2 n_tempCO2 n_tempNit -ascii;
    EESload = ddeexec(tempchan,'[Solve]');

    	T(k) = csvread('tempin.csv');
    ddeterm(tempchan);

    %Calculating new amount of substance
    n_tot(k) = n_tempH2O + n_tempO2 + n_tempCO2 + n_tempNit;
    P(k) = n_tot(k)*R*T(k)/V_chamber;
    
    T_e(k)   = (P_amb./P(k))^(1-1/gamma) * T(k);
    rho_e(k) = P_amb./(R*T_e(k));
    v_e(k)   = sqrt(2*(P(k)-P_amb)./rho_e(k));
    m_out(k) = A_exit*v_e(k)*rho_e(k);
    
    %calculating mass of matter inside rocket
    m_inside(k) = n_tempH2O * M_H2O + n_tempO2 * M_O2 + n_tempCO2 * M_CO2 + n_tempNit * M_Nit;

    %finding change in mass from matter flowing out. m_out has units [kg/s],
    %thus multiplying by dt is required
    delta_m(k)    = m_inside(k) - m_out(k)*dt;

    %Air leaves the rocket first, as it is closest to the nozzle
    n_air = n_air - m_out(k)/M_air * dt;
    
    %---------------------------INITIATE EES TRANSFER--------------------------
    tempchan = ddeinit('EES','DDE');
    %Opens EES work-file
    EESload  = ddeexec(tempchan,'[Open enthalpytotemp2.EES]');

    %Loads in temporary values to save to EES with air abundances reducing per.
    %timestep.
        P_temp = P(k);
        if n_air < 0.05
            n_tempH2O = n_tempH2O;
            n_tempO2  = n_tempO2  - 0.79 * m_out(k)/M_air * dt;
            n_tempCO2 = n_tempCO2;
            n_tempNit = n_tempNit - 0.21 * m_out(k)/M_air * dt;
        else
            n_tempH2O = n_tempH2O - 0.2 * m_out(k)/M_air * dt;
            n_tempO2  = n_tempO2  - 0.2 * m_out(k)/M_air * dt;
            n_tempCO2 = n_tempCO2 - 0.2 * m_out(k)/M_air * dt;
            n_tempNit = n_tempNit - 0.2 * m_out(k)/M_air * dt;
        end
        %H_tempref = n_tempH2O * M_H2O * H_water + n_tempO2 * M_H2O * H_oxygen + n_tempCO2 * H_CO2 %+ n_tempNit * H_water
        H_temp    = H_tot(k-1) + DELTAh_decomposition * n_dot_H2O2_1 * M_H2O2 * dt;
    %Saves variables in file enthalpytotemp.txt in ascii format. 
    save enthalpytotemp.txt P_temp H_temp n_tempH2O n_tempO2 n_tempCO2 n_tempNit -ascii;

    %Loads data into EES and performs the function [solve]
    EESload = ddeexec(tempchan,'[Solve]');
    %EES loads the data, reads vars, calculates T based on enthalpy, saves
    %temperature in tempin.csv which is loaded into T_tot var below.
        T_tot(k) = csvread('tempin.csv');
    %Terminates data exchange
    ddeterm(tempchan);
    H_tot(k) = H_temp;
    P_tot(k) = n_tot(k)*R*T_tot(k)/V_chamber;
    
end

max(H_tot)
max(P_tot)
max(T_tot)
max(n_tot)
n_tempNit
n_tempO2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------Initiate combustion sequence:-----------------------------
j = k
%% 

H_ref_3 = n_dot_H2O_3 * M_H2O * H_water * dt + n_dot_O2_3 * M_O2 * H_oxygen * dt + n_dot_CO2_3 * M_CO2 * H_CO2;
H_3 = H_ref_3 + DELTAh_decomposition * n_dot_H2O2_1 * M_H2O2 * dt + DELTAh_combustion * m_dot_PLA_3 * dt * 15;

%n_tempO2 = 5*n_tempO2
H_tot(k) = H_tot(k) + H_3
for k = j:60
    k = k + 1
    %---------------------------INITIATE EES TRANSFER----------------------
    tempchan = ddeinit('EES','DDE');
%Opens EES work-file
    EESload  = ddeexec(tempchan,'[Open enthalpytotemp2.EES]');

        P_temp = P_tot(k-1);
        H_temp = H_tot(k-1);
        n_tempH2O = n_tempH2O + n_dot_H2O_3 * dt;
        n_tempO2  = n_tempO2 + n_dot_O2_3 * dt;
        n_tempCO2 = n_tempCO2 + n_dot_CO2_3 * dt;
        n_tempNit = n_tempNit;

    save enthalpytotemp.txt P_temp H_temp n_tempH2O n_tempO2 n_tempCO2 n_tempNit -ascii;
    EESload = ddeexec(tempchan,'[Solve]');

    	T(k) = csvread('tempin.csv');
    ddeterm(tempchan);

    %Calculating new amount of substance
    n_tot(k) = n_tempH2O + n_tempO2 + n_tempCO2 + n_tempNit;
    P(k) = n_tot(k)*R*T(k)/V_chamber;
    
    T_e(k)   = (P_amb./P(k))^(1-1/gamma) * T(k);
    rho_e(k) = P_amb./(R*T_e(k));
    v_e(k)   = sqrt(2*(P(k)-P_amb)./rho_e(k));
    m_out(k) = A_exit*v_e(k)*rho_e(k)*100;
    
    %calculating mass of matter inside rocket
    m_inside(k) = n_tempH2O * M_H2O + n_tempO2 * M_O2 + n_tempCO2 * M_CO2 + n_tempNit * M_Nit;

    %finding change in mass from matter flowing out. m_out has units [kg/s],
    %thus multiplying by dt is required
    delta_m(k)    = m_inside(k) - m_out(k)*dt;

    %Air leaves the rocket first, as it is closest to the nozzle
    n_air = n_air - m_out(k)/M_air * dt;
    
    %---------------------------INITIATE EES TRANSFER--------------------------
    tempchan = ddeinit('EES','DDE');
    %Opens EES work-file
    EESload  = ddeexec(tempchan,'[Open enthalpytotemp2.EES]');

    %Loads in temporary values to save to EES with air abundances reducing per.
    %timestep.
%         P_temp = P(k);
%         n_tempH2O = n_tempH2O;
%         n_tempO2  = n_tempO2  - 0.79 * m_out(k)/M_air * dt;
%         n_tempCO2 = n_tempCO2;
%         n_tempNit = n_tempNit - 0.21 * m_out(k)/M_air * dt;
        if n_air < 0.05
            n_tempH2O = n_tempH2O;
            n_tempO2  = n_tempO2  - 0.79 * m_out(k)/M_air * dt;
            n_tempCO2 = n_tempCO2;
            n_tempNit = n_tempNit - 0.21 * m_out(k)/M_air * dt;
        else
            n_tempH2O = n_tempH2O - 0.2 * m_out(k)/M_air * dt;
            n_tempO2  = n_tempO2  - 0.2 * m_out(k)/M_air * dt;
            n_tempCO2 = n_tempCO2 - 0.2 * m_out(k)/M_air * dt;
            n_tempNit = n_tempNit - 0.2 * m_out(k)/M_air * dt;
        end
        %H_ref_3 = n_tempH2O * M_H2O * H_water * dt + n_tempO2 * M_O2 * H_oxygen * dt + n_tempCO2 * M_CO2 * H_CO2;
        H_temp = H_ref_3 + H_tot(k-1) + DELTAh_decomposition * n_dot_H2O2_1 * M_H2O2 * dt + DELTAh_combustion * m_dot_PLA_3 * dt;

    %Saves variables in file enthalpytotemp.txt in ascii format. 
    save enthalpytotemp.txt P_temp H_temp n_tempH2O n_tempO2 n_tempCO2 n_tempNit -ascii;

    %Loads data into EES and performs the function [solve]
    EESload = ddeexec(tempchan,'[Solve]');
    %EES loads the data, reads vars, calculates T based on enthalpy, saves
    %temperature in tempin.csv which is loaded into T_tot var below.
        T_tot(k) = csvread('tempin.csv');
    %Terminates data exchange
    ddeterm(tempchan);
    n_tot(k) = n_tempH2O + n_tempO2 + n_tempCO2 + n_tempNit;
    H_tot(k) = H_temp;
    P_tot(k) = n_tot(k)*R*T_tot(k)/V_chamber;
    
end

%--------------------- Figure Plotting ------------------------------------
total_time = length(P_tot)*dt;
tspan = linspace(0,total_time,length(P_tot));

figure(1)
    yyaxis left
    title('Pressure and temperature over time')
    plot(tspan,P_tot/100,'-o')
    hold on
    axis([0 max(tspan) 1 max(P_tot)/100])
    ylabel('Pressure [bar]')
    yyaxis right 
    plot(tspan,T_tot-T_amb,'-o')
    axis([0 max(tspan) 0 max(T_tot)-273.15])
    ylabel('Temperature [C]')

    xlabel('time [s]')

save('P_totsim.txt', 'P_tot','-ASCII','-append')
save('T_totsim.txt', 'T_tot','-ASCII','-append')

