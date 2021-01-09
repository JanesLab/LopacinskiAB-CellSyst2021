%Figure 5I: Viral antagonism encoded in a lumped-parameter module of type I interferon signaling from autocrine“paracrine sources

clear all
close all hidden

%% Run simulations
rng(1000)
%Simulation with no IFN antiviral response or viral antagonism of the IFN response. Default time changed to 24 hpi for comparison to experimental data.
[NoIFNNoAntMedianSols,NoIFNNoAntMeanSols,NoIFNNoAntUpQuantSols,NoIFNNoAntLowQuantSols,~] = CVB3ODEEval(10,'MaxTime',24,'IFNSwitch','off','VirResponse','off','PlotResults','off');

rng(1000)
%Simulation with just IFN antiviral response but no viral antagonism of the IFN response. Default time changed to 24 hpi for comparison to experimental data
[NoAntMedianSols,NoAntMeanSols,NoAntUpQuantSols,NoAntLowQuantSols,~] = CVB3ODEEval(10,'MaxTime',24,'VirResponse','off','PlotResults','off');

rng(1000)
%Simulation using standard model with default time changed to 24 hpi for comparison to experimental data
[FullModelMedianSols,FullModelMeanSols,FullModelUpQuantSols,FullModelLowQuantSols,~] = CVB3ODEEval(10,'MaxTime',24,'PlotResults','off');

rng(1000)
%Simulation using standard model with IFN co-stimulation at start of infection. Default time changed to 24 hpi for comparison to experimental data
[FullModelCostimMedianSols,FullModelCostimMeanSols,FullModelCostimUpQuantSols,FullModelCostimLowQuantSols,~] = CVB3ODEEval(10,'MaxTime',24,'IFNStimulation','on','IFNStimulationTime',0,'PlotResults','off');

rng(1000)
%Simulation using standard model with IFN stimulation at 6 hpi. Default time changed to 24 hpi for comparison to experimental data
[FullModelStim6MedianSols,FullModelStim6MeanSols,FullModelStim6UpQuantSols,FullModelStim6LowQuantSols,~] = CVB3ODEEval(10,'MaxTime',24,'IFNStimulation','on','IFNStimulationTime',6,'PlotResults','off');

%% Plot Results

%Set up plotting variables for each panel of the subplot
for i = 1:5
    
    if i == 1
        MedianTable = NoIFNNoAntMedianSols;
        UpperQuantileTable = NoIFNNoAntUpQuantSols;
        LowerQuantileTable = NoIFNNoAntLowQuantSols;
    elseif i == 2
        MedianTable = NoAntMedianSols;
        UpperQuantileTable = NoAntUpQuantSols;
        LowerQuantileTable = NoAntLowQuantSols;
    elseif i == 3
        MedianTable = FullModelMedianSols;
        UpperQuantileTable = FullModelUpQuantSols;
        LowerQuantileTable = FullModelLowQuantSols;
    elseif i == 4
        MedianTable = FullModelCostimMedianSols;
        UpperQuantileTable = FullModelCostimUpQuantSols;
        LowerQuantileTable = FullModelCostimLowQuantSols;
    elseif i ==5
        MedianTable = FullModelStim6MedianSols;
        UpperQuantileTable = FullModelStim6UpQuantSols;
        LowerQuantileTable = FullModelStim6LowQuantSols;
    end

%Create the subplots
    Time = table2array(MedianTable(:,{'Time'}));
    
    figure(3)
    subplot(1,5,i)
    Median = table2array(MedianTable(:,{'Virions'}));
    UQ = table2array(UpperQuantileTable(:,{'Virions'}));
    LQ = table2array(LowerQuantileTable(:,{'Virions'}));
    hold on
    patch([Time' fliplr(Time')], [LQ' fliplr(UQ')],'k','EdgeColor','none', 'FaceColor','black','FaceAlpha',.1)
    plot(Time,UQ,'k--','LineWidth',0.5)
    plot(Time,LQ,'k--','LineWidth',0.5)
    plot(Time,Median,'k-','LineWidth',1)
    ylim([0,80])
    xlim([0,25])
    hold off
    
    %Give each subplot a title
    if i == 1
        font_ax(sprintf('%s','No Antiviral Response, No Antagonism'),'Time (h)','Mature Virions (nM)',10,'bold',1.0,0)
    elseif i == 2
        font_ax(sprintf('%s','Antiviral Response, No Antagonism'),'Time (h)','Mature Virions (nM)',10,'bold',1.0,0)
    elseif i == 3
        font_ax(sprintf('%s','Antiviral Response, Antagonism'),'Time (h)','Mature Virions (nM)',10,'bold',1.0,0)
    elseif i == 4
        font_ax(sprintf('%s','Antiviral Response, Antagonism, IFN Stim @ t = 0 hr'),'Time (h)','Mature Virions (nM)',10,'bold',1.0,0)
    elseif i ==5
        font_ax(sprintf('%s','Antiviral Response, Antagonism, IFN Stim @ t = 6 hr'),'Time (h)','Mature Virions (nM)',10,'bold',1.0,0)
    end
  
end

%% Creates function for more concise plot settings

function [] = font_ax(t,xlab,ylab,fsize_ax,fweight,bwidth,cbar_flag)
    %Function: Sets features of the plot display
    %
    %INPUTS:
        %t: Plot title
        %xlab: X-axis label
        %ylab: Y-axis label
        %fsize_ax: Sets the figure's font size
        %fweight: Sets font to bold if desired
        %bwidth: Thickness of the lines in the display
        %cbar_flag: 1 adds a color bar, 0 does not

    title(t,'FontSize',fsize_ax+5,'FontWeight',fweight);
    % title(['Selected Slice (Depth ',num2str(d_ind),'mm)'],'FontSize',15,'FontWeight','bold');
    % xh=xlab;
    % yh=ylab;
    xh=xlabel(xlab);
    yh=ylabel(ylab);
    set([xh,yh],'FontSize',fsize_ax,'FontWeight','bold');
    set(gca,'FontSize',fsize_ax,'FontWeight','bold','LineWidth',bwidth);
    if(cbar_flag==1)
        h=colorbar;
        set(h,'FontSize',fsize_ax,'FontWeight','bold','LineWidth',3);
    end
end