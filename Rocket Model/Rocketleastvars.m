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
%injection duration,
V_chamber   = 0.003;                %[m^3]
P_amb       = 101.3;            %[kPa]
T_amb       = 273.15 + 20;      %[K]
x_oxidizer  = 0.8;              %fraction
f           = 0.59;             %fraction
m_dot_oxidizer = 0.273;         %[kg/s]
mass_H2O2   = 1.3 * x_oxidizer; %[kg]
mass_H2O    = 1.3 - mass_H2O2;  %[kg]
R           = 8.3145;           %[kJ/(mol*K)]
A_chamber   = (0.094)^2;        %[m]^2
A_throat    = (0.0213)^2;       %[m]^2
dt          = 0.01;             %[s]
k           = 2;                %[s]
H_water     = 83.93;            %[kJ/kg]
H_oxygen    = -4.887;           %[kJ/kg]
H_CO2       = -5.168;           %[kJ/kg]
kB          = 1.38065E-23;      %[J/K]
NA          = 6.02214086E23;    %[mol]^-1
time_state1 = (mass_H2O2+mass_H2O)/m_dot_oxidizer; %[s]

%Defining: H2O2, H2O, O2, Plastic grain, CO2 Molar masses:
M_H2O2 = 0.0340147; %[kg/mol]
M_H2O  = 0.0180153; %[kg/mol]
M_O2   = 0.0319988; %[kg/mol]
M_PLA  = 0.0720000; %[kg/mol]
M_CO2  = 0.0440095; %[kg/mol]

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
m_dot_1     = n_dot_H2O_1*M_H2O + n_dot_H2O2_1*M_H2O2;

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
m_dot_2      = n_dot_H2O_2 * M_H2O + n_dot_O2_2 * M_O2; %[kg]

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
m_PLA_3  = n_CO2_3/3*M_PLA;
m_dot_PLA_3  = n_dot_CO2_3/3*M_PLA;
m_3      = n_H2O_3 * M_H2O + n_O2_3 * M_O2 + n_CO2_3 * M_CO2;
m_dot_3  = n_dot_H2O_3 * M_H2O + n_dot_O2_3 * M_O2 + n_dot_CO2_3 * M_CO2;
M_3      = m_3/n_flow_3;



% ---------------------- Total Flow rates ---------------------------------

n_dot_H2O  = [n_dot_H2O_1 n_dot_H2O_2 n_dot_H2O_3];
n_dot_H2O2 = [n_dot_H2O2_1 n_dot_H2O2_2 n_dot_H2O2_3];
n_dot_O2   = [n_dot_O2_1 n_dot_O2_2 n_dot_O2_3];
n_dot_CO2  = [n_dot_CO2_1 n_dot_CO2_2 n_dot_CO2_3];

n_dot_1    = n_dot_H2O2_1 + n_dot_H2O_1 + n_dot_O2_1 + n_dot_CO2_1;
n_dot_2    = n_dot_H2O2_2 + n_dot_H2O_2 + n_dot_O2_2 + n_dot_CO2_2;
n_dot_3    = n_dot_H2O2_3 + n_dot_H2O_3 + n_dot_O2_3 + n_dot_CO2_3;
n_dot_tot  = [n_flow_1 n_flow_2 n_flow_3]*dt;

%Checking that mass is conserved in two first states
mass_conservation = [m_1 m_2 m_3];

%------------Reference enthalpy for further calculations:------------------

H_ref_1 = n_dot_H2O_1 * M_H2O * H_water + n_dot_O2_1 * M_O2 * H_oxygen + n_dot_CO2_1 * M_CO2 * H_CO2
H_ref_2 = n_dot_H2O_2 * M_H2O * H_water + n_dot_O2_2 * M_O2 * H_oxygen + n_dot_CO2_2 * M_CO2 * H_CO2
H_ref_3 = n_dot_H2O_3 * M_H2O * H_water + n_dot_O2_3 * M_O2 * H_oxygen + n_dot_CO2_3 * M_CO2 * H_CO2

H_s2    = H_ref_2 + DELTAh_decomposition*n_dot_H2O2_1*M_H2O2 %[kJ]
H_s3    = H_ref_3 + DELTAh_decomposition*n_dot_H2O2_1*M_H2O2 + DELTAh_combustion*m_dot_PLA_3 %[kJ]

Q_n_dot_3 = H_s3
T_3     = 2/3 * (Q_n_dot_3)/(kB*NA)

P_exit = 92
T_3    = 1920
gamma  = 1.197
M_dot_3  = m_dot_3/n_dot_3
v_exit = sqrt(T_3 * R/M_dot_3 * 2 * gamma/(gamma - 1) * (1-(P_exit/1100)^((gamma-1)/gamma)))
%Thrust calculation
F_exit = m_dot_3*v_exit


% figure(1)
%     title('Pressure per time simulation of injection')
%         
%     subplot(2,2,1);
%         plot(tspan,P_tot,'-o')
%         hold on
%         plot(kbarlinex,kbarliney,'-')
%         hold on
%         plot(tspan,P_outflow,'-o')
%         axis([0 .1 0 10000])
%         xlabel('time [s]')
%         ylabel('Pressure [kPa]')
%         legend('Pressure from initial injection','Working pressure','Pressure with outflow')
%     
%     subplot(2,2,2)
%         plot(tspan,v_out)
%         hold on
%         plot(tspan,max_v)
%         hold on
%         plot(tspan,v_outflow)
%         axis([0 .5 0 4000])
%         xlabel('time [s]')
%         ylabel('Velocity [m/s]')
%         legend('Velocity based on stagnation pressure','Maximum attainable velocity','Velocity with chamber outflow')
%         
% % figure(2)
%      subplot(2,2,3)       
%         %plot(kelvbarlinex,kelvbarliney,'-')
%         hold on
%         plot(T(1:10),P_tot(1:10))
%         hold on
%         axis([0 max(T) 0 max(P_tot)])
%         xlabel('Temperature [K]')
%         ylabel('Pressure [kPa]')
%         legend('Temperature in stagnated chamber')