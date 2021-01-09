%Figure S4A and S4B: Mass balance of capsid, 3Dpol, protease, and 2CATPase under limiting conditions vs. standard model

%% Run simulations

clear all
close all hidden

%Simulation of standard model
rng(1000)
[StandardModelMedianSols,~,StandardModelUpQuantSols,StandardModelLowQuantSols] = CVB3ODEEvalS4Standard(10,'MaxTime',8,'PlotResults','off');

%Simulation with no viral protein or RNA degradation, no defective virions, not counting the capsid protein from infecting virus, and with VirResponse and IFNSwitch set to ‘off’
rng(1000)
[LimConditionsMedianSols,~,LimConditionsUpQuantSols,LimConditionsLowQuantSols] = CVB3ODEEvalS4LimConditions(10,'MaxTime',8,'IFNSwitch','off','VirResponse','off','PlotResults','off');


%% Plot Results

OutputVector =  {'Total Capsid Protein' 'Total 3Dpol' 'Total Protease' 'Total 2C ATPase'};

%Set up plotting variables for each panel of the subplot
for i = 1:2        
    if i == 1
        MedianTable = StandardModelMedianSols;
        UpperQuantileTable = StandardModelUpQuantSols;
        LowerQuantileTable = StandardModelLowQuantSols;
    else
        MedianTable = LimConditionsMedianSols;
        UpperQuantileTable = LimConditionsUpQuantSols;
        LowerQuantileTable = LimConditionsLowQuantSols;
    end
    
    %Create the subplots
    Time = table2array(MedianTable(:,{'Time'}));
    
    subplot(1,2,i)
    for j = 1:4
        Median = table2array(MedianTable(:,{OutputVector{j}}));
        UQ = table2array(UpperQuantileTable(:,{OutputVector{j}}));
        LQ = table2array(LowerQuantileTable(:,{OutputVector{j}}));
        h(j)=semilogy(Time,Median,'r-','LineWidth',2);
        hold on
        switch j
            case 1
                h(j).Color=[195 22 28]/255;
            case 2
                h(j).Color=[247 148 30]/255;
            case 3
                h(j).Color=[3 149 63]/255;
            case 4
                h(j).Color=[144 63 152]/255;
        end
        g=semilogy(Time,LQ,'r--','LineWidth',0.75);
        switch j
            case 1
                g.Color=[195 22 28]/255;
            case 2
                g.Color=[247 148 30]/255;
            case 3
                g.Color=[3 149 63]/255;
            case 4
                g.Color=[144 63 152]/255;
        end
        g=semilogy(Time,UQ,'r--','LineWidth',0.75);
        switch j
            case 1
                g.Color=[195 22 28]/255;
            case 2
                g.Color=[247 148 30]/255;
            case 3
                g.Color=[3 149 63]/255;
            case 4
                g.Color=[144 63 152]/255;
        end
    end
    hold off
    legend(h,{'Total Capsid Protein' 'Total 3Dpol' 'Total Protease' 'Total 2CATPase'})
    
    %Give each subplot a title
    if i == 1
        font_ax(sprintf('%s','Standard Model'),'Time (h)','Concentration (nM)',10,'bold',1.0,0)
    else
        font_ax(sprintf('%s','Limiting Conditions Model'),'Time (h)','Concentration (nM)',10,'bold',1.0,0)
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