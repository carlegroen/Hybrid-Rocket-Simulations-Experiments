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
V_chamber   = 3;                %[L]
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



% ---------------------- Total Flow rates ---------------------------------

n_dot_H2O  = [n_dot_H2O_1 n_dot_H2O_2 n_dot_H2O_3];
n_dot_H2O2 = [n_dot_H2O2_1 n_dot_H2O2_2 n_dot_H2O2_3];
n_dot_O2   = [n_dot_O2_1 n_dot_O2_2 n_dot_O2_3];
n_dot_CO2  = [n_dot_CO2_1 n_dot_CO2_2 n_dot_CO2_3];

n_dot_1    = n_dot_H2O2_1 + n_dot_H2O_1 + n_dot_O2_1 + n_dot_CO2_1;
n_dot_2    = n_dot_H2O2_2 + n_dot_H2O_2 + n_dot_O2_2 + n_dot_CO2_2;
n_dot_3    = n_dot_H2O2_3 + n_dot_H2O_3 + n_dot_O2_3 + n_dot_CO2_3;
n_dot_tot        = [n_flow_1 n_flow_2 n_flow_3]*dt;

%Checking that mass is conserved in two first states
mass_conservation = [m_1 m_2 m_3];

%------------Reference enthalpy for further calculations:------------------
H_ref(1) = n_dot_H2O_1 * M_H2O * H_water + n_dot_O2_1 * M_O2 * H_oxygen + n_dot_CO2_1 * M_CO2 * H_CO2;
H_ref(2) = n_dot_H2O_2 * M_H2O * H_water + n_dot_O2_2 * M_O2 * H_oxygen + n_dot_CO2_2 * M_CO2 * H_CO2;
H_ref(3) = n_dot_H2O_3 * M_H2O * H_water + n_dot_O2_3 * M_O2 * H_oxygen + n_dot_CO2_3 * M_CO2 * H_CO2;

for z=1:3
    H_dot_ref(z) = n_dot_H2O(z) * M_H2O * H_water + n_dot_O2(z) * M_O2 * H_oxygen + n_dot_CO2(z) * M_CO2 * H_CO2;
end
H_comparison = H_dot_ref;

%------------Total temperature and enthalpy for real states:---------------
H_real(1) = H_ref(1);
H_real(2) = H_ref(2) + DELTAh_decomposition*n_dot_H2O2_1*M_H2O2;
H_real(3) = H_ref(3) + DELTAh_decomposition*n_dot_H2O2_1*M_H2O2 + DELTAh_combustion*m_dot_PLA_3;
H_dot_real = H_real*dt;
H_Real_3 = H_real(3);
%Create pressure array
P_tot = [];
T     = [];
H     = [];
v_outflow = [];

%Assign setting at t=dt, just after injection
P_tot(1) = P_amb;
dP = n_dot_1*dt * T_amb * R/V_chamber;
P_tot(2) = P_amb + dP;
T(1)     = T_amb;
T(1) = 275.15
T(2) = 278.15
%assign setting at t=2dt, after first decomposition
dP =(n_dot_2+n_dot_1)*dt * T_amb * R/V_chamber;
%P_tot(3) = P_tot(1) + dP;
%T(2) = T_amb + 2/3 * (H_dot_real(1)*dt+dP*V_chamber)/(kB*n_flow_1*dt*NA);

density   = (mass_H2O2+mass_H2O)/V_chamber;
%v_outflow(2) = sqrt(2*(P_tot(2)-P_amb)/density);
%v_outflow(3) = sqrt(2*(P_tot(3)-P_amb)/density);
%P = P_amb;
%gamma =1.197;

dpP = n_dot_tot(1) * T_amb * R/V_chamber;
T_2 = T_amb + 2/3 * ((H_dot_real(1)+H_dot_real(2))+1100*V_chamber)/(kB*(n_dot_tot(1)+n_dot_tot(2))*NA);
T_2 = 301.1
%%%%%%%Calculating Temperature for state 1%%%%%%%%%%%

chan2   = ddeinit('EES','DDE');
EESload = ddeexec(chan2,'[Open P_CALC2.EES]');
    P1      = 101.3;
    T1      = 293.15;
    save Pressureout.txt P1 T1 H_Real_3 -ascii
    EESload = ddeexec(chan2,'[Solve]');
    H1      = csvread('EESPressureout.csv');
    H_ref1  = n_dot_H2O_1*M_H2O*H1(1) + n_dot_O2_1*M_O2*H1(2) + n_dot_CO2_1*M_CO2*H1(3);
ddeterm(chan2);

chan2   = ddeinit('EES','DDE');
EESload = ddeexec(chan2,'[Open P_CALC2.EES]');
    P1      = 101.3;
    T1      = 293.15;
    save Pressureout.txt P1 T1 H_Real_3 -ascii
    EESload = ddeexec(chan2,'[Solve]');
    H1      = csvread('EESPressureout.csv');
    H_ref2  = n_dot_H2O_2*M_H2O*H1(1)*dt + n_dot_O2_2*M_O2*H1(2)*dt + n_dot_CO2_2*M_CO2*H1(3)*dt + DELTAh_decomposition*n_dot_H2O2_1*M_H2O2*dt;
ddeterm(chan2);
    H(1) = H_ref2;
    H(2) = H_ref2;
    
   size(P_tot)
   size(T)
    
%Initializes data transfer to EES     
chan3   = ddeinit('EES','DDE');
%Opens EES work-file
EESload = ddeexec(chan3,'[Open enthalpytotemp.EES]');  
t= [];
for k = 2:50
         t = k*dt;
         n_out     = [sum(n_dot_H2O) sum(n_dot_O2) sum(n_dot_CO2)]*dt*0.85;
         n_in      = [sum(n_dot_H2O) sum(n_dot_O2) sum(n_dot_CO2)]*dt;
         n_delta   = n_in - n_out;
         
         P_tot(k)  = P_tot(k-1) + sum(n_delta)/V_chamber * R * T(k-1);
         H(k) = H(k-1) + DELTAh_decomposition * n_dot_H2O2_1 * dt * M_H2O2 + DELTAh_combustion * m_dot_PLA_3*dt;
    
         %Loads temporary working pressure and enthalpy
         P_temp = P_tot(k);      
         H_temp = H(k);
    %saves .txt file with data, loads in EES and solves equations.
    save enthalpytotemp.txt P_temp H_temp -ascii;
    EESload = ddeexec(chan3,'[Solve]');
        %Saves EES results in T(k) var
        T(k)    = csvread('tempin.csv');
end

    %closes data transfer
    ddeterm(chan3);

%t=2*dt:dt:time_state1/dt
% t= []
% for k = 2:1000
%          t = k*dt;
%          m_dot_out = 1;
%          n_out     = (n_dot_1+n_dot_2+n_dot_3)*dt*0.75;
%          n_in      = (n_dot_1+n_dot_2+n_dot_3)*dt;
%          n_delta   = n_in - n_out;
%          %k         = k + 1;
%          
%          
%          
%          T(k)      = T(k-1) + 2/3 * ((sum(H_dot_real))*dt+dP*V_chamber)/(kB*n_delta*NA);
%          dP        = n_delta * T(k) * R/V_chamber;
%          P_tot(k)  = P_tot(k-1) + dP;
%          %T(k)      = T(k-1) + 2/3 * ((sum(H_dot_real))*dt+dP*V_chamber)/(kB*(n_dot_1+n_dot_2+n_dot_3)*dt*NA);
%          P_outflow(k) = P_tot(k-1) + n_in * T(k-1) * R/V_chamber;
%          v_outflow(k) = sqrt(2*(P_outflow(k)-P_amb)/density);
% end

v_out = sqrt(2*(P_tot-P_amb)/density);



% ---------------------- Plotting Process ---------------------------------
%Timesteps for process.
tspan = (0:49)*dt/time_state1;

%Lines at specified ranges
kbarlinex = [0 5];
kbarliney = [1000 1000];
kelvbarliney = [1000 1100];
kelvbarlinex = [273.15 2000];

%Approximate velocity-cap
T_test    = 1920;

m_dot_H2O_3 = n_dot_H2O_3 * M_H2O;
m_dot_O2_3  = n_dot_O2_3 * M_O2
m_dot_CO2_3 = n_dot_O2_3 * M_CO2
M_dot_3     = (m_dot_H2O_3 + m_dot_O2_3 + m_dot_CO2_3)/n_dot_3


gamma = 1.2
T_end = 2225+273;
T_sim = linspace(T_amb,T_end,200);
v_max = sqrt(gamma*R*T_sim/M_dot_3);
T_end = 2500;
T_sim = linspace(T_amb,T_end,1000);
P_sim = sqrt(T_sim)*6+900;

figure(1)
    title('Pressure per time simulation of injection')
        
    subplot(2,2,1);
        plot(tspan,P_tot,'-o')
        hold on
        plot(kbarlinex,kbarliney,'-')
        hold on
        %plot(tspan,P_outflow,'-o')
        axis([0 .5 0 5000])
        xlabel('time [s]')
        ylabel('Pressure [kPa]')
        legend('Pressure in chamber','Working pressure','Pressure without mass outflow')
    
    subplot(2,2,2)
        plot((0:49)*dt/time_state1,v_out)
        hold on
        plot((1:200)*dt,v_max)
        hold on
        %plot(tspan,v_outflow)
        axis([0 1 0 1500])
        xlabel('time [s]')
        ylabel('Velocity [m/s]')
        legend('Velocity based on stagnation pressure','Maximum attainable velocity','Velocity with chamber outflow')
        
% figure(2)
     subplot(2,2,3)       
        plot(kelvbarlinex,kelvbarliney,'-')
        hold on
        plot(T(1:36),P_tot(1:36))
        hold on
        plot(T_sim,P_sim)
        hold on
        axis([T_amb 3500 P_amb 1500])
        xlabel('Temperature [K]')
        ylabel('Pressure [kPa]')
        legend('Temperature in stagnated chamber','Simulated Pressure','Theoretical Pressure')