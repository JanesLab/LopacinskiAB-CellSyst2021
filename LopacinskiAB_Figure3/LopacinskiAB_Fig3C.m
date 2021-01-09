clear all
close all hidden

rng(1000); %set pseudorandom number seed

%Load absolute stranded qPCR data
[qPCR,varnames,timepoints]=tblread('CVB3_pos-negstrand_absolutequant.csv', ...
    ',');
timepoints=str2num(timepoints);

rc = 100; %RunCount index to speed troubleshooting
rng(1000)
[med,~,uq,lq]=CVB3ODEEval(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc);
median_pos=med{:,57}; %Total +ssRNA
uq_pos=uq{:,57};
lq_pos=lq{:,57};
median_neg=med{:,58}; %Total -ssRNA
uq_neg=uq{:,58};
lq_neg=lq{:,58};

subplot(1,2,1)
h=semilogy(med.Time,median_pos,'r-','LineWidth',2);
h.Color=[121 49 127]/255;
hold on
h=semilogy(med.Time,uq_pos,'r--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=semilogy(med.Time,lq_pos,'r--','LineWidth',0.75);
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

subplot(1,2,2)
h=semilogy(med.Time,median_neg,'r-','LineWidth',2);
h.Color=[121 49 127]/255;
hold on
h=semilogy(med.Time,uq_neg,'r--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=semilogy(med.Time,lq_neg,'r--','LineWidth',0.75);
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