clear all
close all hidden

rc = 100; %RunCount index to speed troubleshooting

PolysomeRange = linspace(1,30,30);

for i = 1:length(PolysomeRange)
    rng(1000)
    [~,~,~,~,~,VPmin(:,i)]=CVB3ODEEval_FigS3A(10,'MaxTime',8, ...
        'PlotResults','off','RunCount',rc,'PolysomeSize',PolysomeRange(i));
end

h=plot(PolysomeRange, median(VPmin),'r-','LineWidth',2);
h.Color=[35 31 32]/255;
hold on
h=plot(PolysomeRange,quantile(VPmin,0.05),'r--', ...
    'LineWidth',0.75);
h.Color=[209 211 212]/255;
h=plot(PolysomeRange,quantile(VPmin,0.95),'r--', ...
    'LineWidth',0.75);
h.Color=[209 211 212]/255;
xlabel('Number of ribosomes per polysome')
ylabel('Time to exponential protein synthesis (hr)')
title('Effect of polysome size on timing of viral protein synthesis')
axis([0 30 0 10])