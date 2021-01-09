clear all
close all hidden
warning off

rc = 100; %RunCount index to speed troubleshooting

CARvect = linspace(5000,9000,5);
for i=1:length(CARvect)
    rng(1000)
    [med,~,uq,lq]=CVB3ODEEvalFigS2FG(10,'MaxTime',12,'PlotResults','off', ...
        'RunCount',rc,'CAR',CARvect(i),'DAF',4.4e4,'DIP',800); %using DAF estimates from AC16-CAR cells
    median_CAR_VP1(:,i)=med{:,60}; %Total viral protein
    uq_CAR_VP1(:,i)=uq{:,60};
    lq_CAR_VP1(:,i)=lq{:,60};
end

DIPvect = [0 50 100 150 200];
for i=1:length(DIPvect)
    rng(1000)
    [med,~,uq,lq]=CVB3ODEEvalFigS2FG(10,'MaxTime',12,'PlotResults','off', ...
        'RunCount',rc,'CAR',2.3e3,'DAF',6.1e4,'DIP',DIPvect(i)); %using CAR and DAF estimates from AC16 cells
    median_DIP_VP1(:,i)=med{:,60}; %Total viral protein
    uq_DIP_VP1(:,i)=uq{:,60};
    lq_DIP_VP1(:,i)=lq{:,60};
end

subplot(1,2,1)
h=semilogy(med.Time,median_CAR_VP1,'-','LineWidth',2);
h(1).Color=[195 22 28]/255;
h(2).Color=[247 148 30]/255;
h(3).Color=[255 194 14]/255;
h(4).Color=[3 149 63]/255;
h(5).Color=[0 94 158]/255;
hold on
eb_color=reshape([h(1:5).Color],3,5)';
for i=1:length(CARvect)
    semilogy(med.Time,uq_CAR_VP1(:,i),'--','LineWidth',0.75,'Color', ...
        eb_color(i,:));
    semilogy(med.Time,lq_CAR_VP1(:,i),'--','LineWidth',0.75,'Color', ...
        eb_color(i,:));
end
xlabel('Time (hr)')
ylabel('Viral protein (nM)')
axis([0 12 1e-4 1e3])
title('CAR threshold')
legend(num2str(CARvect'))

subplot(1,2,2)
h=semilogy(med.Time,median_DIP_VP1,'-','LineWidth',2);
h(1).Color=[195 22 28]/255;
h(2).Color=[247 148 30]/255;
h(3).Color=[255 194 14]/255;
h(4).Color=[3 149 63]/255;
h(5).Color=[0 94 158]/255;
hold on
eb_color=reshape([h(1:5).Color],3,5)';
for i=1:length(DIPvect)
   semilogy(med.Time,uq_DIP_VP1(:,i),'--','LineWidth',0.75,'Color', ...
       eb_color(i,:));
   semilogy(med.Time,lq_DIP_VP1(:,i),'--','LineWidth',0.75,'Color', ...
       eb_color(i,:));
end
xlabel('Time (hr)')
ylabel('Viral protein (nM)')
axis([0 12 1e-4 1e3])
title('DIP threshold')
legend(num2str(DIPvect'))

