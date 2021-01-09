clear all
close all hidden

rc = 100; %RunCount index to speed troubleshooting

rng(1000)
[med,~,uq,lq]=CVB3ODEEvalFigS2A(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc,'CAR',2.3e3,'Recycling','off');
median_uCAR(:,1,1)=med{:,12}; %Unbound CAR
uq_uCAR(:,1,1)=uq{:,12};
lq_uCAR(:,1,1)=lq{:,12};

rng(1000)
[med,~,uq,lq]=CVB3ODEEvalFigS2A(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc,'CAR',2.3e3,'Recycling','on');
median_uCAR(:,2,1)=med{:,12}; %Unbound CAR
uq_uCAR(:,2,1)=uq{:,12};
lq_uCAR(:,2,1)=lq{:,12};

rng(1000)
[med,~,uq,lq]=CVB3ODEEvalFigS2A(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc,'CAR',5.5e6,'Recycling','off');
median_uCAR(:,1,2)=med{:,12}; %Unbound CAR
uq_uCAR(:,1,2)=uq{:,12};
lq_uCAR(:,1,2)=lq{:,12};

rng(1000)
[med,~,uq,lq]=CVB3ODEEvalFigS2A(10,'MaxTime',8,'PlotResults','off', ...
    'RunCount',rc,'CAR',5.5e6,'Recycling','on');
median_uCAR(:,2,2)=med{:,12}; %Unbound CAR
uq_uCAR(:,2,2)=uq{:,12};
lq_uCAR(:,2,2)=lq{:,12};

subplot(1,2,1)
h=semilogy(med.Time,median_uCAR(:,1,1),'r-','LineWidth',2);
h.Color=[3 149 63]/255;
hold on
h=semilogy(med.Time,median_uCAR(:,2,1),'b-','LineWidth',2);
h.Color=[121 49 127]/255;
h=semilogy(med.Time,uq_uCAR(:,1,1),'r--','LineWidth',0.75);
h.Color=[3 149 63]/255;
h=semilogy(med.Time,lq_uCAR(:,1,1),'r--','LineWidth',0.75);
h.Color=[3 149 63]/255;
h=semilogy(med.Time,uq_uCAR(:,2,1),'b--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=semilogy(med.Time,lq_uCAR(:,2,1),'b--','LineWidth',0.75);
h.Color=[121 49 127]/255;
xlabel('Time (hr)')
ylabel('Unbound CAR (nM)')
title('AC16 cells')
legend({'Degraded','Recycled'})
axis([0 8 1 1e6])

subplot(1,2,2)
h=semilogy(med.Time,median_uCAR(:,1,2),'r-','LineWidth',2);
h.Color=[3 149 63]/255;
hold on
h=semilogy(med.Time,median_uCAR(:,2,2),'b-','LineWidth',2);
h.Color=[121 49 127]/255;
h=semilogy(med.Time,uq_uCAR(:,1,2),'r--','LineWidth',0.75);
h.Color=[3 149 63]/255;
h=semilogy(med.Time,lq_uCAR(:,1,2),'r--','LineWidth',0.75);
h.Color=[3 149 63]/255;
h=semilogy(med.Time,uq_uCAR(:,2,2),'b--','LineWidth',0.75);
h.Color=[121 49 127]/255;
h=semilogy(med.Time,lq_uCAR(:,2,2),'b--','LineWidth',0.75);
h.Color=[121 49 127]/255;
xlabel('Time (hr)')
ylabel('Unbound CAR (nM)')
title('AC16-CAR cells')
legend({'Degraded', 'Recycled'})
axis([0 8 1 1e6])
