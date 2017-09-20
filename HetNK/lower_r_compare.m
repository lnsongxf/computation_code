addpath('export_fig');

% Plot
InCompMkt = load('EqTransAgg.mat');
CompMkt = load('EqTransCompMkt.mat');

figure;clf;hold on;
plot( (log(CompMkt.EqTrans.Y_t(1:40)) - log(CompMkt.EqTrans.Y_t(end))) * 100,'k-','LineWidth',2);
plot( (log(InCompMkt.EqTransAgg.Y_t(1:40)) - log(InCompMkt.EqTransAgg.Y_t(end))) * 100,'r-','LineWidth',2);
legend('Complete Market','Incomplete Market');
xlabel('Quarter');
ylabel('Percentage');
title('Output');
export_fig('graph/lower_r_Y.pdf','-dpdf','-transparent');

figure;clf;hold on;
plot( (CompMkt.EqTrans.Pi_t(1:40)-1) * 100,'k-','LineWidth',2);
plot( (InCompMkt.EqTransAgg.Pi_t(1:40)-1) * 100,'r-','LineWidth',2);
legend('Complete Market','Incomplete Market');
xlabel('Quarter');
ylabel('Percentage');
title('Inflation');
export_fig('graph/lower_r_Pi.pdf','-dpdf','-transparent');