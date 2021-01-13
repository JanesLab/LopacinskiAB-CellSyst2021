% This function regenerates the data for Figures 6A, 6B, S6A, and S6B
% from Lopacinski et al. (2020)
clear all
close all hidden

rc = 100; % sets run count to 100

DetectorDeg = [0.001 0.005];
IFNStimulationTime = [4 5 6];
RunsRemaining = length(DetectorDeg)*(length(IFNStimulationTime)+1);

for i = 1:length(DetectorDeg)
    rng(1000) % set the seed for consistency between runs
    RunsRemaining
    [med,~,uq,lq] = CVB3ODEEval(10,'MaxTime',24,'EC50_DetectorDeg',...
        DetectorDeg(i),'PlotResults','off','RunCount',rc);
    
    median_filledVirions(:,i,1) = med{:,43}; % filled virions
    uq_filledVirions(:,i,1) = uq{:,43};
    lq_filledVirions(:,i,1) = lq{:,43};
    
    median_ISProtein(:,i,1) = med{:,56}; % interferon stimulated proteins
    uq_ISProtein(:,i,1) = uq{:,56};
    lq_ISProtein(:,i,1) = lq{:,56};
    
    median_viralProtease(:,i,1) = med{:,23}; % viral protease
    uq_viralProtease(:,i,1) = uq{:,23};
    lq_viralProtease(:,i,1) = lq{:,23};
    
    RunsRemaining = RunsRemaining - 1;
end

for j = 2:length(IFNStimulationTime)+1   
    for i = 1:length(DetectorDeg)
        rng(1000)
        RunsRemaining
        [med,~,uq,lq] = CVB3ODEEval(10,'MaxTime',24,'EC50_DetectorDeg',...
            DetectorDeg(i),'IFNStimulation','on','IFNStimulationTime',...
            IFNStimulationTime(j-1),'PlotResults','off','RunCount',rc);
        
        median_filledVirions(:,i,j) = med{:,43}; % filled virions
        uq_filledVirions(:,i,j) = uq{:,43};
        lq_filledVirions(:,i,j) = lq{:,43};
        
        median_ISProtein(:,i,j) = med{:,56}; % interferon stimulated proteins
        uq_ISProtein(:,i,j) = uq{:,56};
        lq_ISProtein(:,i,j) = lq{:,56};
        
        median_viralProtease(:,i,j) = med{:,23}; % viral protease
        uq_viralProtease(:,i,j) = uq{:,23};
        lq_viralProtease(:,i,j) = lq{:,23};

        RunsRemaining = RunsRemaining - 1;
    end
end
%%
figure(1) % Figure 6A
subplot(2,2,1)
plot(med.Time,median_filledVirions(:,1,1),'r-','LineWidth',0.5) % DetectorDeg = 0.001, no exogenous IFN
hold on
plot(med.Time,uq_filledVirions(:,1,1),'r--','LineWidth',0.5)
plot(med.Time,lq_filledVirions(:,1,1),'r--','LineWidth',0.5)
plot(med.Time,median_filledVirions(:,2,1),'b-','LineWidth',0.5) % DetectorDeg = 0.005, no exogenous IFN
plot(med.Time,uq_filledVirions(:,2,1),'b--','LineWidth',0.5)
plot(med.Time,lq_filledVirions(:,2,1),'b--','LineWidth',0.5)
axis([0 25 0 100])
xlabel('Time (hr)')
ylabel('Virions (nM)')
title('No paracrine IFN')
legend('Base model','uq','lq','5x resistant')

subplot(2,2,2)
plot(med.Time,median_filledVirions(:,1,2),'r-','LineWidth',0.5) % DetectorDeg = 0.001, IFN @ 4 hr
hold on
plot(med.Time,uq_filledVirions(:,1,2),'r--','LineWidth',0.5)
plot(med.Time,lq_filledVirions(:,1,2),'r--','LineWidth',0.5)
plot(med.Time,median_filledVirions(:,2,2),'b-','LineWidth',0.5) % DetectorDeg = 0.005, IFN @ 4 hr
plot(med.Time,uq_filledVirions(:,2,2),'b--','LineWidth',0.5)
plot(med.Time,lq_filledVirions(:,2,2),'b--','LineWidth',0.5)
axis([0 25 0 100])
xlabel('Time (hr)')
ylabel('Virions (nM)')
title('IFN = 4 hr')

subplot(2,2,3)
plot(med.Time,median_filledVirions(:,1,3),'r-','LineWidth',0.5) % DetectorDeg = 0.001, IFN @ 5 hr
hold on
plot(med.Time,uq_filledVirions(:,1,3),'r--','LineWidth',0.5)
plot(med.Time,lq_filledVirions(:,1,3),'r--','LineWidth',0.5)
plot(med.Time,median_filledVirions(:,2,3),'b-','LineWidth',0.5) % DetectorDeg = 0.005, IFN @ 5hr
plot(med.Time,uq_filledVirions(:,2,3),'b--','LineWidth',0.5)
plot(med.Time,lq_filledVirions(:,2,3),'b--','LineWidth',0.5)
axis([0 25 0 100])
xlabel('Time (hr)')
ylabel('Virions (nM)')
title('IFN = 5 hr')

subplot(2,2,4)
plot(med.Time,median_filledVirions(:,1,4),'r-','LineWidth',0.5) % DetectorDeg = 0.001, IFN @ 6 hr
hold on
plot(med.Time,uq_filledVirions(:,1,4),'r--','LineWidth',0.5)
plot(med.Time,lq_filledVirions(:,1,4),'r--','LineWidth',0.5)
plot(med.Time,median_filledVirions(:,2,4),'b-','LineWidth',0.5) % DetectorDeg = 0.005, IFN @ 6 hr
plot(med.Time,uq_filledVirions(:,2,4),'b--','LineWidth',0.5)
plot(med.Time,lq_filledVirions(:,2,4),'b--','LineWidth',0.5)
axis([0 25 0 100])
xlabel('Time (hr)')
ylabel('Virions (nM)')
title('IFN = 6 hr')

figure(2) % Figure 6B
plot(0.5,median_filledVirions(1440,1,1),'ro','MarkerFaceColor','r') % DetetcDeg = 0.001, no  exogenous IFN
hold on
plot(0.5,uq_filledVirions(1440,1,1),'rx','MarkerSize',14)
plot(0.5,lq_filledVirions(1440,1,1),'rx','MarkerSize',14)
plot(1,median_filledVirions(1440,1,2),'ro','MarkerFaceColor','r') % IFN @ 4
plot(1,uq_filledVirions(1440,1,2),'rx','MarkerSize',14)
plot(1,lq_filledVirions(1440,1,2),'rx','MarkerSize',14)
plot(1.5,median_filledVirions(1440,1,3),'ro','MarkerFaceColor','r') % IFN @ 5
plot(1.5,uq_filledVirions(1440,1,3),'rx','MarkerSize',14)
plot(1.5,lq_filledVirions(1440,1,3),'rx','MarkerSize',14)
plot(2,median_filledVirions(1440,1,4),'ro','MarkerFaceColor','r') % IFN @ 6
plot(2,uq_filledVirions(1440,1,4),'rx','MarkerSize',14)
plot(2,lq_filledVirions(1440,1,4),'rx','MarkerSize',14)

plot(0.5,median_filledVirions(1440,2,1),'bo','MarkerFaceColor','b') % DetectDeg = 0.005, no exogenous IFN
plot(0.5,uq_filledVirions(1440,2,1),'bx','MarkerSize',14)
plot(0.5,lq_filledVirions(1440,2,1),'bx','MarkerSize',14)
plot(1,median_filledVirions(1440,2,2),'bo','MarkerFaceColor','b') % IFN @ 4
plot(1,uq_filledVirions(1440,2,2),'bx','MarkerSize',14)
plot(1,lq_filledVirions(1440,2,2),'bx','MarkerSize',14)
plot(1.5,median_filledVirions(1440,2,3),'bo','MarkerFaceColor','b') % IFN @ 5
plot(1.5,uq_filledVirions(1440,2,3),'bx','MarkerSize',14)
plot(1.5,lq_filledVirions(1440,2,3),'bx','MarkerSize',14)
plot(2,median_filledVirions(1440,2,4),'bo','MarkerFaceColor','b') % IFN @ 6
plot(2,uq_filledVirions(1440,2,4),'bx','MarkerSize',14)
plot(2,lq_filledVirions(1440,2,4),'bx','MarkerSize',14)

axis([0 2.5 0 100])
xlabel('IFN stimulation time (hr)')
ylabel('Virions at 24 hr (nM)')

figure(3) % Figure S6A
subplot(1,2,1)
plot(0.5,median_ISProtein(1440,1,1),'ro','MarkerFaceColor','r') % DetectDeg = 0.001, no  exogenous IFN
hold on
plot(0.5,uq_ISProtein(1440,1,1),'rx','MarkerSize',14)
plot(0.5,lq_ISProtein(1440,1,1),'rx','MarkerSize',14)
plot(0.5,median_ISProtein(1440,2,1),'bo','MarkerFaceColor','b') % DetectDeg = 0.005, no exogenous IFN
plot(0.5,uq_ISProtein(1440,2,1),'bx','MarkerSize',14)
plot(0.5,lq_ISProtein(1440,2,1),'bx','MarkerSize',14)

axis([0 1 0 350])
ylabel('IFN Stimulated Protein (nM)')
xlabel('No supplemental IFN')
legend('base','uq','lq','5x resistant')

subplot(1,2,2)
plot(1,median_ISProtein(1440,1,2),'ro','MarkerFaceColor','r') % DetectDeg = 0.001, IFN @ 4
hold on
plot(1,uq_ISProtein(1440,1,2),'rx','MarkerSize',14)
plot(1,lq_ISProtein(1440,1,2),'rx','MarkerSize',14)
plot(1.5,median_ISProtein(1440,1,3),'ro','MarkerFaceColor','r') % IFN @ 5
plot(1.5,uq_ISProtein(1440,1,3),'rx','MarkerSize',14)
plot(1.5,lq_ISProtein(1440,1,3),'rx','MarkerSize',14)
plot(2,median_ISProtein(1440,1,4),'ro','MarkerFaceColor','r') % IFN @ 6
plot(2,uq_ISProtein(1440,1,4),'rx','MarkerSize',14)
plot(2,lq_ISProtein(1440,1,4),'rx','MarkerSize',14)
plot(1,median_ISProtein(1440,2,2),'bo','MarkerFaceColor','b') % DetectDeg = 0.005, IFN @ 4
plot(1,uq_ISProtein(1440,2,2),'bx','MarkerSize',14)
plot(1,lq_ISProtein(1440,2,2),'bx','MarkerSize',14)
plot(1.5,median_ISProtein(1440,2,3),'bo','MarkerFaceColor','b') % IFN @ 5
plot(1.5,uq_ISProtein(1440,2,3),'bx','MarkerSize',14)
plot(1.5,lq_ISProtein(1440,2,3),'bx','MarkerSize',14)
plot(2,median_ISProtein(1440,2,4),'bo','MarkerFaceColor','b') % IFN @ 6
plot(2,uq_ISProtein(1440,2,4),'bx','MarkerSize',14)
plot(2,lq_ISProtein(1440,2,4),'bx','MarkerSize',14)

axis([0.5 2.5 6000 19000])
ylabel('IFN Stimulated Protein (nM)')
xlabel('Supplemental IFN')

figure(4) % Figure S6B
subplot(1,2,1)
plot(0.5,median_viralProtease(1440,1,1),'ro','MarkerFaceColor','r') % DetectDeg = 0.001, no  exogenous IFN
hold on
plot(0.5,uq_viralProtease(1440,1,1),'rx','MarkerSize',14)
plot(0.5,lq_viralProtease(1440,1,1),'rx','MarkerSize',14)
plot(0.5,median_viralProtease(1440,2,1),'bo','MarkerFaceColor','b') % DetectDeg = 0.005, no exogenous IFN
plot(0.5,uq_viralProtease(1440,2,1),'bx','MarkerSize',14)
plot(0.5,lq_viralProtease(1440,2,1),'bx','MarkerSize',14)

axis([0 1 0 2000])
ylabel('Viral Protease (nM)')
xlabel('No supplemental IFN')
legend('base','uq','lq','5x resistant')

subplot(1,2,2)
plot(1,median_viralProtease(1440,1,2),'ro','MarkerFaceColor','r') % DetectDeg = 0.001, IFN @ 4
hold on
plot(1,uq_viralProtease(1440,1,2),'rx','MarkerSize',14)
plot(1,lq_viralProtease(1440,1,2),'rx','MarkerSize',14)
plot(1.5,median_viralProtease(1440,1,3),'ro','MarkerFaceColor','r') % IFN @ 5
plot(1.5,uq_viralProtease(1440,1,3),'rx','MarkerSize',14)
plot(1.5,lq_viralProtease(1440,1,3),'rx','MarkerSize',14)
plot(2,median_viralProtease(1440,1,4),'ro','MarkerFaceColor','r') % IFN @ 6
plot(2,uq_viralProtease(1440,1,4),'rx','MarkerSize',14)
plot(2,lq_viralProtease(1440,1,4),'rx','MarkerSize',14)
plot(1,median_viralProtease(1440,2,2),'bo','MarkerFaceColor','b') % DetectDeg = 0.005, IFN @ 4
plot(1,uq_viralProtease(1440,2,2),'bx','MarkerSize',14)
plot(1,lq_viralProtease(1440,2,2),'bx','MarkerSize',14)
plot(1.5,median_viralProtease(1440,2,3),'bo','MarkerFaceColor','b') % IFN @ 5
plot(1.5,uq_viralProtease(1440,2,3),'bx','MarkerSize',14)
plot(1.5,lq_viralProtease(1440,2,3),'bx','MarkerSize',14)
plot(2,median_viralProtease(1440,2,4),'bo','MarkerFaceColor','b') % IFN @ 6
plot(2,uq_viralProtease(1440,2,4),'bx','MarkerSize',14)
plot(2,lq_viralProtease(1440,2,4),'bx','MarkerSize',14)

axis([0.5 2.5 0 20])
ylabel('Viral Protease (nM)')
xlabel('Supplemental IFN')
