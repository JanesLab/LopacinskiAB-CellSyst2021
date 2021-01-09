%% Figure S4C and S4D.  Population-level calibration of VP1 and eIF4G cleavage

clear all
close all hidden
rng(1000)

%Import experimental data. AC16-CAR cells were infected at an MOI of 10 and immunoblotted for VP1 expression and eIF4G cleavage at 0, 1, 2, 4, 5, 6, and 8 hours post-infection.
experimentalData = readtable('CVB3_proteins_relativequant.csv');

%Replace non-detected values with 0. Assume the species was not present if not detectable.
experimentalDataVars = experimentalData{:,2:end};
experimentalDataVars(isnan(experimentalDataVars))=0;
experimentalData{:,2:end}=experimentalDataVars;

% Execute simulations
[~,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',8,'PlotResults','off','PopulationSetting','Population');

%Estimate the adjustment constant to convert pixel intensity to approximate intracellular concentration. 
%Adjustment constant will be estimated using linear regression. Regressionmodel : meanModelPrediction ~ adjustmentCoefficient * experimentalData

%Find the mean scaling coefficient for VP1 timecourse data
modelTimePointIndexes = [1,61,121,241,301,361,480];
VP1DataScaling = experimentalData{:,2:5}./meanResultsTable{modelTimePointIndexes,60};
adjustmentCoefficientVP1 = mean(VP1DataScaling,'all','omitnan');

%Find a global adjustment parameter for eIF4G timecourse data
adjFitClveIF4G = fitlm(reshape(experimentalData{:,10:13},height(experimentalData(:,10:13))*width(experimentalData(:,10:13)),1), ...
    repmat(meanResultsTable{modelTimePointIndexes,'Viral Ribosomes'} - meanResultsTable{1,'Viral Ribosomes'}, ...
    width(experimentalData(:,10:13)),1),'Intercept',false);
adjustmentCoefficientClveIF4G=adjFitClveIF4G.Coefficients{1,1};
    
%Adjust the western blot data from pixels to VP1 and protease intracellular concentration.
adjustedVP1Data = experimentalData{:,2:5}./adjustmentCoefficientVP1;
adjustedClveIF4G = meanResultsTable{1,'Viral Ribosomes'} + experimentalData{:,10:13}*adjustmentCoefficientClveIF4G;

%Find the geometric mean of the experimental data
adjVP1DataMean = geomean(adjustedVP1Data,2);
adjClveIF4GDataMean = geomean(adjustedClveIF4G,2);

%Find the log-transformed standard error
adjVP1DataSE = std(adjustedVP1Data')'/sqrt(length(adjustedVP1Data));
adjClveIF4GDataSE = std(adjustedClveIF4G')'/sqrt(length(adjustedClveIF4G));

% Plot simulations overlayed with adjusted experimental results
subplot(1,2,1)
 modelTime = meanResultsTable{:,1};
 modelmean = meanResultsTable{:,60};
 modelUQ = Q95ResultsTable{:,60};
 modelLQ = Q05ResultsTable{:,60};
 
 %remove absolute zeros from model and data so that it can plotted on a semilog axis.
 modelmean(modelmean==0)=1e-4;
 modelUQ(modelUQ==0)=1e-4;
 modelLQ(modelLQ==0)=1e-4;
 
 semilogy(modelTime,modelmean,'Color','#7E2F8E','LineWidth',1)
 ylim([0.0001,10000])
 hold on
 patch([modelTime' fliplr(modelTime')], [modelLQ' fliplr(modelUQ')],'k','EdgeColor','none', 'FaceColor','#7E2F8E','FaceAlpha',.25)
 semilogy(experimentalData{:,1},adjustedVP1Data,'k.','MarkerSize',15)
 errorbar(experimentalData{:,1},adjVP1DataMean,adjVP1DataSE,'k.','LineWidth',1)
 semilogy(experimentalData{:,1},adjVP1DataMean,'s','MarkerEdgeColor','none','MarkerFaceColor','#77AC30','MarkerSize',5)
 font_ax(sprintf('%s','Viral Protein VP1'),'Time (h)','Concentration (nM)',10,'bold',1.0,0)
 hold off
 
 subplot(1,2,2)
 modelTime = meanResultsTable{:,1};
 modelmean = meanResultsTable{:,'Viral Ribosomes'};
 modelUQ = Q95ResultsTable{:,'Viral Ribosomes'};
 modelLQ = Q05ResultsTable{:,'Viral Ribosomes'};
 
 %remove absolute zeros from model and data so that it can plotted on a semilog axis.
 modelmean(modelmean==0)=1e-6;
 modelUQ(modelUQ==0)=1e-6;
 modelLQ(modelLQ==0)=1e-6;
 
 plot(modelTime,modelmean,'Color','#7E2F8E','LineWidth',1)
 ylim([4,11])
 hold on
 patch([modelTime' fliplr(modelTime')], [modelLQ' fliplr(modelUQ')],'k','EdgeColor','none', 'FaceColor','#7E2F8E','FaceAlpha',.25)
 plot(experimentalData{:,1},adjustedClveIF4G,'k.','MarkerSize',15)
 errorbar(experimentalData{:,1},adjClveIF4GDataMean,adjClveIF4GDataSE,'k.','LineWidth',1)
 plot(experimentalData{:,1},adjClveIF4GDataMean,'s','MarkerEdgeColor','none','MarkerFaceColor','#77AC30','MarkerSize',5)
 font_ax(sprintf('%s','Accessible Ribosomes'),'Time (h)','Concentration (nM)',10,'bold',1.0,0)
 hold off
 
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