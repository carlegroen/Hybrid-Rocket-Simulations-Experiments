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
dt          = 0.01;             %[s]
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
m_PLA_3  = n_CO2_3/3*M_PLA;
m_dot_PLA_3  = n_dot_CO2_3/3*M_PLA;
m_3      = n_H2O_3 * M_H2O + n_O2_3 * M_O2 + n_CO2_3 * M_CO2;
M_3      = m_3/n_flow_3;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------


%Isentropic relations and identities are used to calculate pressure.

dx  = 1/dt;
n_1 = [n_dot_H2O_1; n_dot_O2_1; n_dot_CO2_1]/dx;
n_2 = [n_dot_H2O_2; n_dot_O2_2; n_dot_CO2_2]/dx;
n_3 = [n_dot_H2O_3; n_dot_O2_3; n_dot_CO2_3]/dx;

m_tot1   = n_dot_H2O2_1/dx*M_H2O2 + n_1(1)* M_H2O + n_1(2) * M_O2 + n_1(3)*M_CO2 ;
m_tot2   = n_2(1)* M_H2O + n_2(2) * M_O2 + n_2(3)*M_CO2 ;
m_tot3   = n_3(1)* M_H2O + n_3(2) * M_O2 + n_3(3)*M_CO2% + m_dot_PLA_3/dx;
m_in     = m_tot3;

%mass fractions of exhaust:
m_frac1 = [0.8; (n_1(1)* M_H2O)/m_tot1; (n_1(2)* M_O2)/m_tot1; (n_1(3)* M_CO2)/m_tot1];
m_frac3 = [(n_3(1)* M_H2O)/m_tot3; (n_3(2)* M_O2)/m_tot3; (n_3(3)* M_CO2)/m_tot3; (m_dot_PLA_3/dx)/m_tot3];


%Reference enthalpies

H_ref_2 = n_2(1) * M_H2O * H_water + n_2(2) * M_O2 * H_oxygen + n_2(3) * M_CO2 * H_CO2;
H_ref_3 = n_3(1) * M_H2O * H_water + n_3(2) * M_O2 * H_oxygen + n_3(3) * M_CO2 * H_CO2;

H_2 = H_ref_2 + DELTAh_decomposition * n_dot_H2O2_1/dx * M_H2O2;
H_3 = H_ref_3 + DELTAh_decomposition * n_dot_H2O2_1/dx * M_H2O2 + DELTAh_combustion*m_dot_PLA_3/dx;


%--------------------- t = 0 ---------------------------------------------
%Initial state calculation of temperature and pressure:
%First, temperature and pressure at t = 0 is T_amb and P_amb
T = [];
T(1) = T_amb;
P(1) = P_amb;

%number of substance is calculated as:
n_start = P(1)*V_chamber/(T(1)*R);
n_tot(1)= n_start + sum(n_1)+n_dot_H2O2_1/dx;

T_tot(1) = T(1);
P_tot(1) = n_tot(1)*R*T_tot(1)/V_chamber;



%--------------------- t = 1*dt -------------------------------------------

%Matter has now decomposed, still no combustion:
 
%---------------------------INITIATE EES TRANSFER--------------------------
tempchan = ddeinit('EES','DDE');
%Opens EES work-file
EESload  = ddeexec(tempchan,'[Open enthalpytotemp2.EES]');

    P_temp = P_tot(1);
    H_temp = H_2;
    n_tempH2O = n_2(1);
    n_tempO2  = n_2(2);
    n_tempCO2 = n_2(3);
    n_tempAir = n_start;

save enthalpytotemp.txt P_temp H_temp n_tempH2O n_tempO2 n_tempCO2 n_tempAir -ascii;
EESload = ddeexec(tempchan,'[Solve]');

    T(2) = csvread('tempin.csv');
ddeterm(tempchan);
% %---------------------------KILL EES TRANSFER------------------------------
n_tot(2) = n_start + sum(n_2);
P(2) = n_tot(2)*R*T(2)/V_chamber;
H(2) = H_temp;

T_e(2)   = (P_amb./P(2))^(1-1/gamma) * T(2);
rho_e(2) = P_amb./(R*T_e(2));
v_e(2)   = sqrt(2*(P(2)-P_amb)./rho_e(2));
m_out(2) = A_exit*v_e(2)*rho_e(2)/dx;

delta_m(2) = m_tot2-m_out(2);

%residual air to materials:
n_air_frac = n_start/n_tot(2);
n_mat_frac = sum(n_2)/n_tot(2);

m_frac2 = [(n_2(1)* M_H2O)/m_tot2; (n_2(2)* M_O2)/m_tot2; (n_2(3)* M_CO2)/m_tot2];

tempchan = ddeinit('EES','DDE');
EESload  = ddeexec(tempchan,'[Open enthalpytotemp3.EES]');


    m_H2O_s2 = delta_m(2)*m_frac2(1);
    m_O2_s2  = delta_m(2)*m_frac2(2);
    m_CO2_s2 = delta_m(2)*m_frac2(3);
    m_air_s2 = delta_m(2)*n_air_frac;

    P_temp = P(2);
    H_tempref = m_H2O_s2 * H_water + m_O2_s2 * H_oxygen + m_CO2_s2 * H_CO2
    H_temp = H_tempref + DELTAh_decomposition * n_dot_H2O2_1/dx * M_H2O2

save enthalpytotemp.txt P_temp H_temp m_H2O_s2 m_O2_s2 m_CO2_s2 m_air_s2 -ascii;
EESload = ddeexec(tempchan,'[Solve]');

    T_tot(2) = csvread('tempin.csv');
    
ddeterm(tempchan);
    H_tot(2) = H_ref_3;
    m_CO2_s2 = n_3(3)*M_CO2
    P_tot(2) = n_tot(2)*R*T_tot(2)/V_chamber;


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------- t = 2*dt -------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

T_new=[];

%Fresh attempt!


H_tot(2) = H_3
for k=3:15
tempchan3 = ddeinit('EES','DDE');
%Opens EES work-file
EESload  = ddeexec(tempchan3,'[Open enthalpytotemp2.EES]');
    P_temp = P_tot(k-1);
    H_temp = H_tot(k-1);
    n_tempH2O = n_3(1);
    n_tempO2  = n_3(2);
    n_tempCO2 = n_3(3);
    n_tempAir = n_tempAir*n_air_frac

save enthalpytotemp.txt P_temp H_temp n_tempH2O n_tempO2 n_tempCO2 n_tempAir -ascii;
EESload = ddeexec(tempchan3,'[Solve]');

    T(k) = csvread('tempin.csv');
    
ddeterm(tempchan3);
% %---------------------------KILL EES TRANSFER------------------------------
n_tot(k) = n_tempAir + sum(n_3);
P(k) = n_tot(k)*R*T(k)/V_chamber;
H(k) = H_temp;

T_e(k)   = (P_amb./P(k))^(1-1/gamma) * T(k);
rho_e(k) = P_amb./(R*T_e(k))*100;
v_e(k)   = sqrt(2*(P(k)-P_amb)./rho_e(k));
m_out(k) = A_exit*v_e(k)*rho_e(k)/dx;

m_tot3   = n_3(1)* M_H2O + n_3(2) * M_O2 + n_3(3)*M_CO2;

delta_m(k) = m_tot3-m_out(k);

%residual air to materials:
n_air_frac = n_tempAir/n_tot(k);
n_mat_frac = sum(n_3)/n_tot(k);

tempchan4 = ddeinit('EES','DDE');
EESload  = ddeexec(tempchan4,'[Open enthalpytotemp3.EES]');


    m_H2O_s2 = delta_m(k)*m_frac3(1);
    m_O2_s2  = delta_m(k)*m_frac3(2);
    m_CO2_s2 = delta_m(k)*m_frac3(3);
    m_air_s2 = delta_m(k)*n_air_frac;

    P_temp = P(k);
    H_tempref = m_H2O_s2 * H_water + m_O2_s2 * H_oxygen + m_CO2_s2 * H_CO2;
    H_temp = H_tempref + DELTAh_decomposition * n_dot_H2O2_1/dx * M_H2O2 + DELTAh_combustion * (m_CO2_s2/M_CO2)/3*M_PLA;

save enthalpytotemp.txt P_temp H_temp m_H2O_s2 m_O2_s2 m_CO2_s2 m_air_s2 -ascii;
EESload = ddeexec(tempchan4,'[Solve]');

    T_tot(k) = csvread('tempin.csv');
    
ddeterm(tempchan4);
    H_tot(k) = H_temp;
    P_tot(k) = n_tot(k)*R*T_tot(k)/V_chamber;
    
    n_3(1) = m_H2O_s2/M_H2O;
    n_3(2)  = m_O2_s2/M_O2;
    n_3(3) = m_CO2_s2/M_CO2;
    n_tempAir = m_air_s2/M_air;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tempchan = ddeinit('EES','DDE');
% Opens EES work-file
% EESload  = ddeexec(tempchan,'[Open enthalpytotemp5.EES]');
% 
% for k = 3:3
%     P_temp = P_tot(k-1);
%     H_temp = H_tot(k-1) + DELTAh_decomposition * n_dot_H2O2_1/dx * M_H2O2 + DELTAh_combustion * (m_CO2_s2/M_CO2)/3*M_PLA;
%     n_tempH2O = n_3(1);
%     n_tempO2  = n_3(2);
%     n_tempCO2 = n_3(3);
%     n_tempAir = m_air_s2/(M_air*k^(1.5));
%     m_PLA     = (m_CO2_s2/M_CO2)/3*M_PLA;
% save enthalpytotemp.txt P_temp H_temp n_tempH2O n_tempO2 n_tempCO2 n_tempAir m_PLA -ascii;
% EESload = ddeexec(tempchan,'[Solve]');
% 
%     T(k) = csvread('tempin.csv');
%   
% tempchan = ddeinit('EES','DDE');
% 
% n_tot(k) = n_start+sum(n_3)+m_CO2_s2/M_CO2;
% P(k) = n_tot(k)*R*T(k)/V_chamber;
% %%%%%%From the other direction%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% T_e(k)   = (P_amb./P(k))^(1-1/gamma) * T(k);
% rho_e(k) = P_amb./(R*T_e(k));
% v_e(k)   = sqrt(2*(P(k)-P_amb)./rho_e(k));
% m_out(k) = A_exit*v_e(k)*rho_e(k)/dx;
% 
% %%%%%%%%Now both sides
% 
% Opens EES work-file
% tempchan2 = ddeinit('EES','DDE');
% EESload  = ddeexec(tempchan2,'[Open enthalpytotemp4.EES]');
% 
%     delta_m(k) = (m_in-m_out(k));
%     
%     n_air_frac = n_tempAir/(n_tot(k)*k^2);
%     n_mat_frac = 1-n_air_frac;
%     m_H2O_s2 = delta_m(k)*m_frac3(1)*n_mat_frac;
%     m_O2_s2  = delta_m(k)*m_frac3(2)*n_mat_frac;
%     m_CO2_s2 = delta_m(k)*m_frac3(3)*n_mat_frac;
%     m_air_s2 = delta_m(k)*n_air_frac;
%     
%     
%     n_tot(k) = n_start + m_H2O_s2/M_H2O + m_O2_s2/M_O2 + m_CO2_s2/M_CO2 + m_air_s2/M_air;
%     
%     P_temp = P(k);
%     
%     H_temp2(2) = m_H2O_s2 * H_water + m_O2_s2 * H_oxygen + m_CO2_s2 * H_CO2;
%     H_temp2(k) = H_temp2(k-1) + DELTAh_decomposition * n_dot_H2O2_1/dx * M_H2O2 + DELTAh_combustion * (m_CO2_s2/M_CO2)/3*M_PLA;
% 
% save enthalpytotemp.txt P_temp H_temp2 m_H2O_s2 m_O2_s2 m_CO2_s2 m_air_s2 -ascii;
% EESload = ddeexec(tempchan,'[Solve]');
% 
%     T_tot(k) = csvread('tempin.csv');
%     
% tempchan2 = ddeinit('EES','DDE');
% 
%     P_tot(k) = n_tot(k)*R*T(k)/V_chamber;
%     H_tot(k) = H_temp2(k);
%     n_3(1) = m_H2O_s2/M_H2O;
%     n_3(2) = m_O2_s2/M_O2;
%     n_3(3) = m_CO2_s2/M_CO2;
%     m_CO2_s2 = m_CO2_s2/M_CO2;
%     m_air_s2 = m_air_s2;
%     
% end
% ddeterm(tempchan2);
% ddeterm(tempchan);

%--------------------- Figure Plotting ------------------------------------

m_dot_H2O_3 = n_dot_H2O_3 * M_H2O;
m_dot_O2_3  = n_dot_O2_3 * M_O2
m_dot_CO2_3 = n_dot_O2_3 * M_CO2
M_dot_3     = (m_dot_H2O_3 + m_dot_O2_3 + m_dot_CO2_3)/sum(n_3)


gamma = 1.2
T_end = 2225+273;
T_sim = linspace(T_amb,T_end,length(P_tot));
v_ = sqrt(gamma*R*T_sim/M_dot_3);
P_sim = sum(n_3)*R.*T_sim./V_chamber 

tspan = linspace(1*dt,length(P_tot)*dt,length(P_tot));


figure(1)
    title('Pressure per time simulation of injection')
        
    subplot(2,1,1);
    	yyaxis left
        plot(tspan,P_tot/100,'-o')
        hold on
        plot(tspan,P_sim/100,'-ok')
        hold on
        yyaxis left
        axis([0 max(tspan) 0 max(P_tot)/100])
        ylabel('Pressure [bar]')
        
        yyaxis right
        plot(tspan,T_tot-T_amb)
        axis([0 max(tspan) 0 max(T_tot)])
        ylabel('Temperature [C]')
        
        xlabel('time [s]')
        legend('Pressure in chamber','Pressure without mass outflow')
    subplot(2,1,2);
    	yyaxis left
        plot(T_tot,P_tot/100,'-o')
        hold on
%     subplot(2,2,2)
%         plot((0:49)*dt/time_state1,v_out)
%         hold on
%         plot((1:200)*dt,v_max)
%         hold on
%         %plot(tspan,v_outflow)
%         axis([0 1 0 1500])
%         xlabel('time [s]')
%         ylabel('Velocity [m/s]')
%         legend('Velocity based on stagnation pressure','Maximum attainable velocity','Velocity with chamber outflow')

