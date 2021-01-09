%Figure S5F: Parameter sensitivity analysis for mature virions at 24 hr and 10 PFU

%% Run simulations and plot results

clear all
close all hidden
warning off
rng(1000)

%Simulation of standard model at 24 hpi with 10 pfu. Sensitivity analysis on and showing virions with default scaling.
[~,~,~,~,SensitivitySolutions] = CVB3ODEEval(10,'MaxTime',24,'PlotResults','off', ...
    'SensitivityAnalysis','on','SensitivityAnalysisOutput','Virion','ScalingFactor',2);

caxis([-3 3])
subplot(1,10,1:9)
caxis([-3 3])