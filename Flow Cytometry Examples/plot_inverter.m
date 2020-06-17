load('analysis/TetR-batch.mat');

Ara_uM = [1000 100 50 10 1 0.5 0.1 0];
ZeroOnLog_Ara_uM = [1000 100 50 10 1 0.5 0.1 0.01];
n_conditions = numel(Ara_uM);

% GFP input, RFP output
GFP_means = nan(n_conditions,1);
GFP_stds = GFP_means; RFP_means = GFP_means; RFP_stds = GFP_means; 
for i=1:n_conditions
    GFP_means(i) = results{i}.means(1);
    GFP_stds(i) = results{i}.stds(1);
    RFP_means(i) = results{i}.means(2);
    RFP_stds(i) = results{i}.stds(2);
end


% Plot induction curve
h = figure('PaperPosition',[1 1 5 5]);
loglog(ZeroOnLog_Ara_uM, GFP_means,'-'); hold on;
LU = geo_error_bars(GFP_means',GFP_stds');
errorbar(ZeroOnLog_Ara_uM, GFP_means, LU(1,:), LU(2,:));
ZeroOnLog(ZeroOnLog_Ara_uM(end),0.03);
xlim([0.007, 2e3]); ylim([3e2 2e5]);
xlabel('L-arabinose (uM)');
ylabel('GFP (MEFL)');
title('pBAD induction');
outputfig(h,'pBAD_induction','plots');

% Plot inverter curve
h = figure('PaperPosition',[1 1 5 5]);
loglog(GFP_means, RFP_means,'-'); hold on;
LU = geo_error_bars(RFP_means',RFP_stds');
errorbar(GFP_means, RFP_means, LU(1,:), LU(2,:));
xlabel('Input (MEFL)');
xlim([1e3 1e5]); ylim([2e1 3e3]);
ylabel('Output (MEFL)');
title('TetR inverter');
outputfig(h,'TetR_inverter','plots');
