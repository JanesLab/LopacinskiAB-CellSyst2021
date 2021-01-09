clear all
close all hidden

rc = 100; %RunCount index to speed troubleshooting

rng(1000)
[med,~,uq,lq]=CVB3ODEEval(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc);
median_VP1(:,1)=med{:,60}; %Total viral protein
uq_VP1(:,1)=uq{:,60};
lq_VP1(:,1)=lq{:,60};

rng(1000)
[med,~,uq,lq]=CVB3ODEEval_noRNAdeg(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc);
median_VP1(:,2)=med{:,60}; %Total viral protein
uq_VP1(:,2)=uq{:,60};
lq_VP1(:,2)=lq{:,60};

rng(1000)
[med,~,uq,lq]=CVB3ODEEval_nodsRNAsens(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc);
median_VP1(:,3)=med{:,60}; %Total viral protein
uq_VP1(:,3)=uq{:,60};
lq_VP1(:,3)=lq{:,60};

rng(1000)
[med,~,uq,lq]=CVB3ODEEval_VRO250(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc);
median_VP1(:,4)=med{:,60}; %Total viral protein
uq_VP1(:,4)=uq{:,60};
lq_VP1(:,4)=lq{:,60};

rng(1000)
[med,~,uq,lq]=CVB3ODEEval_VRO100(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc);
median_VP1(:,5)=med{:,60}; %Total viral protein
uq_VP1(:,5)=uq{:,60};
lq_VP1(:,5)=lq{:,60};

for i=1:4
    subplot(1,4,i)
    h=plot(med.Time,median_VP1(:,i),'r-','LineWidth',2);
    h.Color=[35 31 32]/255;
    hold on
    h=plot(med.Time,uq_VP1(:,i),'r--','LineWidth',0.75);
    h.Color=[209 211 212]/255;
    h=plot(med.Time,lq_VP1(:,i),'r--','LineWidth',0.75);
    h.Color=[209 211 212]/255;
    if i==4
        h=plot(med.Time,median_VP1(:,i+1),'b-','LineWidth',2);
        h.Color=[35 31 32]/255;
        h=plot(med.Time,uq_VP1(:,i+1),'b--','LineWidth',0.75);
        h.Color=[209 211 212]/255;
        h=plot(med.Time,lq_VP1(:,i+1),'b--','LineWidth',0.75);
        h.Color=[209 211 212]/255;
    end
    switch i
        case 1
            title('Control')
        case 2
            title('No RNA degradation')
        case 3
            title('No dsRNA sensing')
        case 4
            title('Reduced VRO effect')
    end
    xlabel('Time (hr)')
    ylabel('Viral protein (nM)')
    axis([0 8 0 300])
end
