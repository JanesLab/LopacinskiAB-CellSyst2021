clear all
close all hidden
warning off

rc = 100; %RunCount index to speed troubleshooting

Efficiency = [0 0.1 0.3 1];
rng(1000)
[Time,PosNegRatio(:,:,1)] = CVB3ODEEvalFig3E(10,'MaxTime',8, ...
    'PlotResults','off','RunCount',100);
h=semilogy(Time,PosNegRatio(:,1,1),'r-','LineWidth',2);
h.Color=[195 22 28]/255;
hold on
h=semilogy(Time,PosNegRatio(:,2,1),'r--','LineWidth',0.5);
h.Color=[195 22 28]/255;
h=semilogy(Time,PosNegRatio(:,3,1),'r--','LineWidth',0.5);
h.Color=[195 22 28]/255;

for i=2:length(Efficiency)
    rng(1000)
    [Time,PosNegRatio(:,:,i)]=CVB3ODEEvalFig3E(10,'MaxTime',8, ...
        'PlotResults','off','RunCount',100,'dsRNAMeasureEfficiency', ...
        Efficiency(i));
    h=semilogy(Time,PosNegRatio(:,1,i),'-','LineWidth',2,'Color', ...
        [0 Efficiency(i) 0.5]);
    switch i
        case 2
            h.Color=[247 148 30]/255;
        case 3
            h.Color=[255 194 14]/255;
        case 4
            h.Color=[3 149 63]/255;
    end
    hold on
    h=semilogy(Time,PosNegRatio(:,2,i),'--','LineWidth',0.5,'Color', ...
        [0 Efficiency(i) 0.5]);
    switch i
        case 2
            h.Color=[247 148 30]/255;
        case 3
            h.Color=[255 194 14]/255;
        case 4
            h.Color=[3 149 63]/255;
    end
    h=semilogy(Time,PosNegRatio(:,3,i),'--','LineWidth',0.5,'Color', ...
        [0 Efficiency(i) 0.5]);
    switch i
        case 2
            h.Color=[247 148 30]/255;
        case 3
            h.Color=[255 194 14]/255;
        case 4
            h.Color=[3 149 63]/255;
    end
end

axis([2 8 1 1e4])
xlabel('Time (hr)')
ylabel('Positive-negative strand ratio')
legend({'0% efficiency' '0% CI (90%)' '0% CI (90%)' ...
    '10% efficiency' '10% CI (90%)' '10% CI (90%)' ...
    '30% efficiency' '30% CI (90%)' '30% CI (90%)' ...
    '100% efficiency' '100% CI (90%)' '100% CI (90%)'})