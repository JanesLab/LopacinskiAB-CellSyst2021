clear all
close all hidden

rc = 100; %RunCount index to speed troubleshooting

rng(1000)
[med,~,uq,lq]=CVB3ODEEvalFig2GH(10,'MaxTime',12,'PlotResults','off', ...
    'RunCount',rc,'CAR',5.5e6,'DAF',5.25e4,'DIP',800);
median_VP1(:,1)=med{:,60}; %Total viral protein
uq_VP1(:,1)=uq{:,60};
lq_VP1(:,1)=lq{:,60};

rng(1000)
[med,~,uq,lq]=CVB3ODEEvalFig2GH(10,'MaxTime',12,'PlotResults','off', ...
    'RunCount',rc,'CAR',2.3e3,'DAF',5.25e4,'DIP',800);
median_VP1(:,2)=med{:,60}; %Total viral protein
uq_VP1(:,2)=uq{:,60};
lq_VP1(:,2)=lq{:,60};

rng(1000)
[med,~,uq,lq]=CVB3ODEEvalFig2GH(10,'MaxTime',12,'PlotResults','off', ...
    'RunCount',rc,'CAR',5.5e6,'DAF',5.25e4,'DIP',0);
median_VP1(:,3)=med{:,60}; %Total viral protein
uq_VP1(:,3)=uq{:,60};
lq_VP1(:,3)=lq{:,60};

rng(1000)
[med,~,uq,lq]=CVB3ODEEvalFig2GH(10,'MaxTime',12,'PlotResults','off', ...
    'RunCount',rc,'CAR',2.3e3,'DAF',5.25e4,'DIP',0);
median_VP1(:,4)=med{:,60}; %Total viral protein
uq_VP1(:,4)=uq{:,60};
lq_VP1(:,4)=lq{:,60};

for i=1:4
    subplot(1,4,i)
    h=semilogy(med.Time,median_VP1(:,i),'r-','LineWidth',2);
    h.Color=[35 31 32]/255;
    hold on
    h=semilogy(med.Time,uq_VP1(:,i),'r--','LineWidth',0.75);
    h.Color=[209 211 212]/255;
    h=semilogy(med.Time,lq_VP1(:,i),'r--','LineWidth',0.75);
    h.Color=[209 211 212]/255;
    xlabel('Time (hr)')
    ylabel('Viral protein (nM)')
    switch i
        case 1
            title('AC16-CAR cells')
        case 2
            title('AC16 cells')
        case 3
            title('AC16-CAR cells, no DIPs')
        case 4
            title('AC16 cells, no DIPs')
    end
    axis([0 12 1e-4 1e4])
end