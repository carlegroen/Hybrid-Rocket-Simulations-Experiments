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
 tspan2 = linspace(2.28,2.48,length(P_tot));
     
 figure(2)
    plot(tspan,B3PEPF*3)
    hold on
    plot(tspan2,P_tot/100,'LineWidth',5)
    hold on
    axis([2.25 3 0 15])
%figure(3)
    %plot(T_tot,P_tot)