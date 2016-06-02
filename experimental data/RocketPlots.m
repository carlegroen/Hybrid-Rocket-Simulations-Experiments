Burn1 = load('Burn1_crop.mat');

    B1FL = UntitledForceLink.Data;

    B1PEPB = UntitledPEpressureback.Data;
    B1PEPI = UntitledPEpressureinjector.Data;
    B1PEPF = UntitledPEpressurefront.Data;

    B1DPB  = UntitledDanfossback.Data;
    B1DPI1 = UntitledDanfossinject1.Data;
    B1DPI2 = UntitledDanfossinject2.Data;
    B1DPF  = UntitledDanfossfront.Data;

    B1PRPF = UntitledPRpressurefront.Data;

    B1Valve = UntitledValvecontrol.Data;
%%
Burn2 = load('Burn2_crop.mat','UntitledForceLink');

    B2FL = UntitledForceLink.Data;

    B2PEPB = UntitledPEpressureback.Data;
    B2PEPI = UntitledPEpressureinjector.Data;
    B2PEPF = UntitledPEpressurefront.Data;

    B2DPB  = UntitledDanfossback.Data;
    B2DPI1 = UntitledDanfossinject1.Data;
    B2DPI2 = UntitledDanfossinject2.Data;
    B2DPF  = UntitledDanfossfront.Data;

    B2PRPF = UntitledPRpressurefront.Data;

    B2Valve = UntitledValvecontrol.Data;
%%
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
    
%% 
    
Burn4 = load('Burn4_crop.mat');


    B4FL = UntitledForceLink.Data;

    B4PEPB = UntitledPEpressureback.Data;
    B4PEPI = UntitledPEpressureinjector.Data;
    B4PEPF = UntitledPEpressurefront.Data;

    B4DPB  = UntitledDanfossback.Data;
    B4DPI1 = UntitledDanfossinject1.Data;
    B4DPI2 = UntitledDanfossinject2.Data;
    B4DPF  = UntitledDanfossfront.Data;

    B4PRPF = UntitledPRpressurefront.Data;

    B4Valve = UntitledValvecontrol.Data;

tspan = linspace(0,10,length(UntitledPEpressureback.Data));
%% 

figure(1)
    plot(tspan,B1PEPB)
    hold on
    plot(tspan,B1PEPI)
    hold on
    plot(tspan,B1PEPF)
    hold on
    plot(tspan,B1DPB)
    hold on
    plot(tspan,B1DPI1)
    hold on
    plot(tspan,B1DPI2)
    hold on
    plot(tspan,B1DPF)
    hold on
    plot(tspan,B1PRPF)
    hold on
  %  plot(tspan,B1Valve)
    
figure(2)
    plot(tspan,B2PEPB)
    hold on
    plot(tspan,B2PEPI)
    hold on
    plot(tspan,B2PEPF)
    hold on
    plot(tspan,B2DPB)
    hold on
    plot(tspan,B2DPI1)
    hold on
    plot(tspan,B2DPI2)
    hold on
    plot(tspan,B2DPF)
    hold on
    plot(tspan,B2PRPF)
    hold on
   % plot(tspan,B2Valve)
    
figure(3)
    plot(tspan,B3PEPB)
    hold on
    plot(tspan,B3PEPI)
    hold on
    plot(tspan,B3PEPF)
    hold on
    plot(tspan,B3DPB)
    hold on
    plot(tspan,B3DPI1)
    hold on
    plot(tspan,B3DPI2)
    hold on
    plot(tspan,B3DPF)
    hold on
    plot(tspan,B3PRPF)
    hold on
   % plot(tspan,B3Valve)
    
figure(4)
    plot(tspan,B4PEPB)
    hold on
    plot(tspan,B4PEPI)
    hold on
    plot(tspan,B4PEPF)
    hold on
    plot(tspan,B4DPB)
    hold on
    plot(tspan,B4DPI1)
    hold on
    plot(tspan,B4DPI2)
    hold on
    plot(tspan,B4DPF)
    hold on
    plot(tspan,B4PRPF*2)
    hold on
    %plot(tspan,B4Valve)
    
    
    figure(5)
    plot(tspan,B1FL,'k')
    hold on
    plot(tspan,B2FL,'y')
    hold on
    plot(tspan,B3FL,'b')
    hold on
    plot(tspan,B4FL,'r')
    hold on