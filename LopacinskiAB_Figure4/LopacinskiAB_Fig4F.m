%% Figure 4F:  Sensitivity to Pentamer:VRO affinity
clear all
close all hidden
rng(1000)

set(groot, 'defaultAxesTickLabelInterpreter','none')
set(groot, 'defaultColorbarTickLabelInterpreter','none')
set(groot, 'defaultGraphplotInterpreter','none')
set(groot, 'defaultLegendInterpreter','none')
set(groot, 'defaultTextInterpreter','none')

%Simulate an MOI of 10 infection for 24hrs.
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ]=CVB3ODEEval(10,'MaxTime',24,'PlotResults','off');

%Replace any zeros with 1e-6 so it can be plotted in the log scale correctly
medianResultsTableVars=medianResultsTable{:,2:end};
Q95ResultsTableVars=Q95ResultsTable{:,2:end};
Q05ResultsTableVars=Q05ResultsTable{:,2:end};

medianResultsTableVars(medianResultsTableVars==0)=1e-6;
Q95ResultsTableVars(Q95ResultsTableVars==0)=1e-6;
Q05ResultsTableVars(Q05ResultsTableVars==0)=1e-6;

medianResultsTable{:,2:end} = medianResultsTableVars;
Q95ResultsTable{:,2:end}=Q95ResultsTableVars;
Q05ResultsTable{:,2:end}=Q05ResultsTableVars;

%Simulate an MOI of 10 infection with the pentamer exiting the VRO increased 10 fold.
rng(1000)
[medianResultsTable_4F1,meanResultsTable_4F1,Q95ResultsTable_4F1,Q05ResultsTable_4F1, ~ ]=CVB3ODEEval_Fig4F1(10,'MaxTime',24,'PlotResults','off');

%Replace any zeros with 1e-6 so it can be plotted in the log scale correctly
medianResultsTable_4F1Vars=medianResultsTable_4F1{:,2:end};
Q95ResultsTable_4F1Vars=Q95ResultsTable_4F1{:,2:end};
Q05ResultsTable_4F1Vars=Q05ResultsTable_4F1{:,2:end};

medianResultsTable_4F1Vars(medianResultsTable_4F1Vars==0)=1e-6;
Q95ResultsTable_4F1Vars(Q95ResultsTable_4F1Vars==0)=1e-6;
Q05ResultsTable_4F1Vars(Q05ResultsTable_4F1Vars==0)=1e-6;

medianResultsTable_4F1{:,2:end} = medianResultsTable_4F1Vars;
Q95ResultsTable_4F1{:,2:end}=Q95ResultsTable_4F1Vars;
Q05ResultsTable_4F1{:,2:end}=Q05ResultsTable_4F1Vars;

%Simulate an MOI of 10 infection with the pentamer exiting the VRO decreased 10 fold.
rng(1000)
[medianResultsTable_4F2,meanResultsTable_4F2,Q95ResultsTable_4F2,Q05ResultsTable_4F2, ~ ]=CVB3ODEEval_Fig4F2(10,'MaxTime',24,'PlotResults','off');

%Replace any zeros with 1e-6 so it can be plotted in the log scale correctly
medianResultsTable_4F2Vars=medianResultsTable_4F2{:,2:end};
Q95ResultsTable_4F2Vars=Q95ResultsTable_4F2{:,2:end};
Q05ResultsTable_4F2Vars=Q05ResultsTable_4F2{:,2:end};

medianResultsTable_4F2Vars(medianResultsTable_4F2Vars==0)=1e-6;
Q95ResultsTable_4F2Vars(Q95ResultsTable_4F2Vars==0)=1e-6;
Q05ResultsTable_4F2Vars(Q05ResultsTable_4F2Vars==0)=1e-6;

medianResultsTable_4F2{:,2:end} = medianResultsTable_4F2Vars;
Q95ResultsTable_4F2{:,2:end}=Q95ResultsTable_4F2Vars;
Q05ResultsTable_4F2{:,2:end}=Q05ResultsTable_4F2Vars;

% Plot simulation results
figure(1)

%Virions
plot(medianResultsTable{:,1},medianResultsTable{:,'Virions'},'Color','#C4151C','LineWidth',1)
xlim([0,24])
ylim([0,80])
hold on
plot(medianResultsTable_4F1{:,1},medianResultsTable_4F1{:,'Virions'},'Color','#FEC10D','LineWidth',1)
plot(medianResultsTable_4F2{:,1},medianResultsTable_4F2{:,'Virions'},'Color','#0089CF','LineWidth',1)
patch([medianResultsTable{:,1}' fliplr(medianResultsTable{:,1}')], [Q05ResultsTable{:,'Virions'}' fliplr(Q95ResultsTable{:,'Virions'}')],'k','EdgeColor','none', 'FaceColor','#C4151C','FaceAlpha',.25)
patch([medianResultsTable_4F1{:,1}' fliplr(medianResultsTable_4F1{:,1}')], [Q05ResultsTable_4F1{:,'Virions'}' fliplr(Q95ResultsTable_4F1{:,'Virions'}')],'k','EdgeColor','none', 'FaceColor','#FEC10D','FaceAlpha',.25)
patch([medianResultsTable_4F2{:,1}' fliplr(medianResultsTable_4F2{:,1}')], [Q05ResultsTable_4F2{:,'Virions'}' fliplr(Q95ResultsTable_4F2{:,'Virions'}')],'k','EdgeColor','none', 'FaceColor','#0089CF','FaceAlpha',.25)
hold off
legend('Complete Model','10x k_Pentamer_off','0.1x k_Pentamer_off')
font_ax(sprintf('%s','Mature Virions'),'Time (h)','Concentration (nM)',10,'bold',1.0,0)

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