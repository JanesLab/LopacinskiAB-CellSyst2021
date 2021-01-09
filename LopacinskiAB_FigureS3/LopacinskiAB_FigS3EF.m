clear all
close all hidden

rng(1000); %set pseudorandom number seed

%Load absolute stranded qPCR data
[qPCR,varnames,timepoints]=tblread('CVB3_pos-negstrand_absolutequant.csv', ...
    ',');
timepoints=str2num(timepoints);

rc = 100; %RunCount index to speed troubleshooting
rng(1000)
[DefTranslate_med,~,DefTranslate_uq,DefTranslate_lq]=CVB3ODEEvalDefTranslate(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc);
DefTranslate_median_pos=DefTranslate_med{:,58}; %Total +ssRNA
DefTranslate_uq_pos=DefTranslate_uq{:,58};
DefTranslate_lq_pos=DefTranslate_lq{:,58};
DefTranslate_median_neg=DefTranslate_med{:,59}; %Total -ssRNA
DefTranslate_uq_neg=DefTranslate_uq{:,59};
DefTranslate_lq_neg=DefTranslate_lq{:,59};

figure(1)
subplot(1,2,1)
h=semilogy(DefTranslate_med.Time,DefTranslate_median_pos,'r-','LineWidth',2);
h.Color=[121 49 127]/255;
hold on
h=semilogy(DefTranslate_med.Time,DefTranslate_uq_pos,'r--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=semilogy(DefTranslate_med.Time,DefTranslate_lq_pos,'r--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=errorbar(timepoints,geomean(qPCR(:,1:4),2),exp(mean(log(qPCR(:,1:4)),2)+ ...
    std(log(qPCR(:,1:4)),0,2)/sqrt(size(qPCR(:,1:4),2)))-geomean(qPCR(:,1:4),2), ...
    geomean(qPCR(:,1:4),2)-exp(mean(log(qPCR(:,1:4)),2)-std(log(qPCR(:,1:4)),0,2)/ ...
    sqrt(size(qPCR(:,1:4),2))),'b+','LineWidth',2);
h.Color=[3 149 63]/255;
h=plot(timepoints,qPCR(:,1:4),'ko');
for i=1:4
    h(i).Color=[35 31 32]/255;
    h(i).MarkerFaceColor=[35 31 32]/255;
end
axis([0 8 1e-3 1e4])
xlabel('Time (hr)')
ylabel('Positive strands (nM)')
title('Defective genomes form')

subplot(1,2,2)
h=semilogy(DefTranslate_med.Time,DefTranslate_median_neg,'r-','LineWidth',2);
h.Color=[121 49 127]/255;
hold on
h=semilogy(DefTranslate_med.Time,DefTranslate_uq_neg,'r--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=semilogy(DefTranslate_med.Time,DefTranslate_lq_neg,'r--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=errorbar(timepoints,geomean(qPCR(:,5:8),2),exp(mean(log(qPCR(:,5:8)),2)+ ...
    std(log(qPCR(:,5:8)),0,2)/sqrt(size(qPCR(:,5:8),2)))-geomean(qPCR(:,5:8),2), ...
    geomean(qPCR(:,5:8),2)-exp(mean(log(qPCR(:,5:8)),2)-std(log(qPCR(:,5:8)),0,2)/ ...
    sqrt(size(qPCR(:,5:8),2))),'b+','LineWidth',2);
h.Color=[3 149 63]/255;
h=plot(timepoints,qPCR(:,5:8),'ko');
for i=1:4
    h(i).Color=[35 31 32]/255;
    h(i).MarkerFaceColor=[35 31 32]/255;
end
axis([0 8 1e-3 1e4])
xlabel('Time (hr)')
ylabel('Negative strands (nM)')
title(' translation complexes')

rng(1000); %set pseudorandom number seed

%Load absolute stranded qPCR data
[qPCR,varnames,timepoints]=tblread('CVB3_pos-negstrand_absolutequant.csv', ...
    ',');
timepoints=str2num(timepoints);

rc = 100; %RunCount index to speed troubleshooting
rng(1000)
[DefTranslateRNALoss_med,~,DefTranslateRNALoss_uq,DefTranslateRNALoss_lq]=CVB3ODEEvalDefTranslateRNALoss(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc);
DefTranslateRNALoss_median_pos=DefTranslateRNALoss_med{:,58}; %Total +ssRNA
DefTranslateRNALoss_uq_pos=DefTranslateRNALoss_uq{:,58};
DefTranslateRNALoss_lq_pos=DefTranslateRNALoss_lq{:,58};
DefTranslateRNALoss_median_neg=DefTranslateRNALoss_med{:,59}; %Total -ssRNA
DefTranslateRNALoss_uq_neg=DefTranslateRNALoss_uq{:,59};
DefTranslateRNALoss_lq_neg=DefTranslateRNALoss_lq{:,59};

figure(2)
subplot(1,2,1)
h=semilogy(DefTranslateRNALoss_med.Time,DefTranslateRNALoss_median_pos,'r-','LineWidth',2);
h.Color=[121 49 127]/255;
hold on
h=semilogy(DefTranslateRNALoss_med.Time,DefTranslateRNALoss_uq_pos,'r--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=semilogy(DefTranslateRNALoss_med.Time,DefTranslateRNALoss_lq_pos,'r--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=errorbar(timepoints,geomean(qPCR(:,1:4),2),exp(mean(log(qPCR(:,1:4)),2)+ ...
    std(log(qPCR(:,1:4)),0,2)/sqrt(size(qPCR(:,1:4),2)))-geomean(qPCR(:,1:4),2), ...
    geomean(qPCR(:,1:4),2)-exp(mean(log(qPCR(:,1:4)),2)-std(log(qPCR(:,1:4)),0,2)/ ...
    sqrt(size(qPCR(:,1:4),2))),'b+','LineWidth',2);
h.Color=[3 149 63]/255;
h=plot(timepoints,qPCR(:,1:4),'ko');
for i=1:4
    h(i).Color=[35 31 32]/255;
    h(i).MarkerFaceColor=[35 31 32]/255;
end
axis([0 8 1e-3 1e4])
xlabel('Time (hr)')
ylabel('Positive strands (nM)')
title('Defective genomes degraded')

subplot(1,2,2)
h=semilogy(DefTranslateRNALoss_med.Time,DefTranslateRNALoss_median_neg,'r-','LineWidth',2);
h.Color=[121 49 127]/255;
hold on
h=semilogy(DefTranslateRNALoss_med.Time,DefTranslateRNALoss_uq_neg,'r--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=semilogy(DefTranslateRNALoss_med.Time,DefTranslateRNALoss_lq_neg,'r--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=errorbar(timepoints,geomean(qPCR(:,5:8),2),exp(mean(log(qPCR(:,5:8)),2)+ ...
    std(log(qPCR(:,5:8)),0,2)/sqrt(size(qPCR(:,5:8),2)))-geomean(qPCR(:,5:8),2), ...
    geomean(qPCR(:,5:8),2)-exp(mean(log(qPCR(:,5:8)),2)-std(log(qPCR(:,5:8)),0,2)/ ...
    sqrt(size(qPCR(:,5:8),2))),'b+','LineWidth',2);
h.Color=[3 149 63]/255;
h=plot(timepoints,qPCR(:,5:8),'ko');
for i=1:4
    h(i).Color=[35 31 32]/255;
    h(i).MarkerFaceColor=[35 31 32]/255;
end
axis([0 8 1e-3 1e4])
xlabel('Time (hr)')
ylabel('Negative strands (nM)')
title(' by nonsense-mediated decay')