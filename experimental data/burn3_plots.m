close all

tspan = linspace(0,10,length(UntitledPEpressureback.Data));
Burn3 = load('Burn3_crop.mat');

    B3FL = UntitledForceLink.Data;

    B3PEPB = UntitledPEpressureback.Data;
    B3PEPI = UntitledPEpressureinjector.Data;
    B3PEPF = UntitledPEpressurefront.Data;

    B3DPB  = UntitledDanfossback.Data;
    B3DPI1 = UntitledDanfossinject1.Data;
    B3DPI2 = UntitledDanfossinject2.Data;
    B3DPF  = UntitledDanfossfront.Data;

    B3PRPF = UntitledPRpressurefront.Data;

    B3Valve = UntitledValvecontrol.Data;
      %plot(tspan,B3FL*2)
    %hold on
    %plot(tspan,B3PEPB)
    %hold on  
figure(1)
     plot(tspan,B3PEPI*4)
    hold on
     plot(tspan,B3PEPF*4)
     hold on
     plot(tspan,B3DPB)
     hold on
     plot(tspan,B3DPI1)
     hold on
     plot(tspan,B3DPI2)
     hold on
     plot(tspan,B3DPF)
     hold on
     plot(tspan,B3PRPF*5)
     hold on
     plot(tspan,B3Valve-5)
     legend('B3PEPI','B3PEPF','B3DPB','B3DPI1','B3DPI2','B3DPF','B3PRPF','B3Valve')
     
     
 T_tot = importdata('T_totsim.txt');
 P_tot = importdata('P_totsim.txt');
 tspan2 = linspace(2.29,2.49,length(P_tot));
  
 
n_3 =          0.1517830457071962
V_chamber = 3;
R           = 8.3145;           %[kJ/(mol*K)]
T_amb = 273.15+20;    
gamma = 1.2;
T_end = 2225+273;
T_sim = linspace(T_amb,T_end,200);
%v_max = sqrt(gamma*R*T_sim/M_dot_3);
T_end = 3500;
T_sim = linspace(T_amb,T_end,1000);
P_sim = n_3*R*T_sim/V_chamber;
 
 
 
 figure(2)
    plot(tspan,B3PEPF*3)
    hold on
    plot(tspan2,P_tot/100,'LineWidth',5)
    hold on
    axis([2.25 3 0 15])
    xlabel('Time [s]')
    ylabel('Pressure [bar]')
    legend('Piezoelectric Pressure measurement of front chamber','Simulated chamber pressure')
figure(3)
    plot(T_tot,P_tot/100)
    xlabel('Temperature [K]')
    ylabel('Pressure [bar]')
    legend('Pressure per temperature')
    
T_tot20 = importdata('T_totsim20.txt');
T_tot60 = importdata('T_totsim60.txt');
T_tot100 = importdata('T_totsim100.txt');
T_tot140 = importdata('T_totsim140.txt');
T_tot180 = importdata('T_totsim180.txt');
T_tot220 = importdata('T_totsim220.txt');
T_tot260 = importdata('T_totsim260.txt');
T_tot300 = importdata('T_totsim300.txt');


P_tot20 = importdata('P_totsim20.txt')/100;
P_tot60 = importdata('P_totsim60.txt')/100;
P_tot100 = importdata('P_totsim100.txt')/100;
P_tot140 = importdata('P_totsim140.txt')/100;
P_tot180 = importdata('P_totsim180.txt')/100;
P_tot220 = importdata('P_totsim220.txt')/100;
P_tot260 = importdata('P_totsim260.txt')/100;
P_tot300 = importdata('P_totsim300.txt')/100;

vexitstd = importdata('vexitstd.txt')

figure(4)
    plot(tspan,B3PEPF*3)
    hold on
        plot(tspan2,P_tot300,'LineWidth',2)
    hold on
        plot(tspan2,P_tot260,'LineWidth',2)
    hold on
        plot(tspan2,P_tot220,'LineWidth',2)
    hold on
        plot(tspan2,P_tot180,'LineWidth',2)
    hold on
        plot(tspan2,P_tot140,'LineWidth',2)
    hold on
        plot(tspan2,P_tot100,'LineWidth',2)
    hold on
         plot(tspan2,P_tot60,'LineWidth',2)
    hold on
        plot(tspan2,P_tot20,'LineWidth',2)
    hold on
        plot(tspan2,5*vexitstd/350,'LineWidth',2)






    axis([2.25 3 0 20])
    xlabel('Time [s]')
    ylabel('Pressure [bar]')
    legend('Measured Pressure','T_{auto} = 300^o C','T_{auto} = 260^o C','T_{auto} = 220^o C','T_{auto} = 180^o C','T_{auto} = 140^o C','T_{auto} = 100^o C','T_{auto} = 60^o C','T_{auto} = 20^o C')
    
 figure(5)
    plot(T_tot20,P_tot20)
    hold on
    plot(T_tot60,P_tot60)
    hold on
    plot(T_tot100,P_tot100)
    hold on
    plot(T_tot140,P_tot140)
    hold on
    plot(T_tot180,P_tot180)
    hold on
    plot(T_tot220,P_tot220)
    hold on
    plot(T_tot260,P_tot260)
    hold on
    plot(T_tot300,P_tot300)
    hold on
    plot(T_sim,P_sim/100)
    xlabel('Temperature [K]')
    ylabel('Pressure [bar]')
    legend('T_{auto} = 20^o C','T_{auto} = 60^o C','T_{auto} = 100^o C','T_{auto} = 140^o C','T_{auto} = 180^o C','T_{auto} = 220^o C','T_{auto} = 260^o C','T_{auto} = 300^o C')   

P_dot01 = importdata('P_totsimmdot0.1.txt')/100;
P_dot02 = importdata('P_totsimmdot0.2.txt')/100;
P_dot03 = importdata('P_totsimmdot0.3.txt')/100;
P_dot04 = importdata('P_totsimmdot0.4.txt')/100;
P_dot05 = importdata('P_totsimmdot0.5.txt')/100;
P_dot015 = importdata('P_totsimmdot0.15.txt')/100;
P_dot025 = importdata('P_totsimmdot0.25.txt')/100;
P_dot035 = importdata('P_totsimmdot0.35.txt')/100;
P_dot045 = importdata('P_totsimmdot0.45.txt')/100;

figure(6)
    plot(tspan,B3PEPF*3)
    hold on
    plot(tspan2,P_dot01,'LineWidth',2)
    hold on
    plot(tspan2,P_dot015,'LineWidth',2)
    hold on
    plot(tspan2,P_dot02,'LineWidth',2)
    hold on
    plot(tspan2,P_dot025,'LineWidth',2)
    hold on
    plot(tspan2,P_dot03,'LineWidth',2)
    hold on
    plot(tspan2,P_dot035,'LineWidth',2)
    hold on
    plot(tspan2,P_dot04,'LineWidth',2)
    hold on
    plot(tspan2,P_dot045,'LineWidth',2)
    hold on
    plot(tspan2,P_dot05,'LineWidth',2)
    hold on
    axis([2.25 3 0 20])
    
    xlabel('Time [s]')
    ylabel('Pressure [bar]')
    legend('Measured P/t','mass flow = 0.1 kg s^{-1}','mass flow = 0.15 kg s^{-1}','mass flow = 0.2 kg s^{-1}','mass flow = 0.25 kg s^{-1}','mass flow = 0.3 kg s^{-1}','mass flow = 0.35 kg s^{-1}','mass flow = 0.4 kg s^{-1}','mass flow = 0.45 kg s^{-1}','mass flow = 0.5 kg s^{-1}')
