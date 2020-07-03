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
xlim([0.007, 2e3]); ylim([3e2 2e5]);
ZeroOnLog(ZeroOnLog_Ara_uM(end),0.03);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute and plot SNR
h = figure('PaperPosition',[1 1 6 4]);
semilogx(results{1}.bincenters,results{1}.bincounts(:,1)./max(results{1}.bincounts(:,1)),'r-'); hold on;
semilogx(results{end}.bincenters,results{end}.bincounts(:,1)./max(results{end}.bincounts(:,1)),'b-');
autofluorescence_boundary = 2*getStdERF(get_autofluorescence_model(CM,1));
plot([autofluorescence_boundary autofluorescence_boundary],[0 2],'k--');
legend('Ara High','Ara Low','AF Limit');
xlim([1e2 1e6]); ylim([0 1.1]);
xlabel('GFP (MEFL)'); ylabel('Normalized Density');
title('pBAD induction plus/minus');
outputfig(h','pBAD_induction_plusminus','plots');

pBAD_SNR = SNR(results{1}.means(1),results{end}.means(1),results{1}.stds(1),results{end}.stds(1))
%    5.7252 dB

h = figure('PaperPosition',[1 1 6 4]);
semilogx(results{1}.bincenters,results{1}.bincounts(:,2)./max(results{1}.bincounts(:,2)),'r-'); hold on;
semilogx(results{end}.bincenters,results{end}.bincounts(:,2)./max(results{end}.bincounts(:,2)),'b-');
autofluorescence_boundary = 2*getStdERF(get_autofluorescence_model(CM,2));
plot([autofluorescence_boundary autofluorescence_boundary],[0 2],'k--');
legend('Location','NorthWest','TetR High','TetR Low','AF Limit');
xlim([1e0 1e4]); ylim([0 1.1]);
xlabel('RFP (MEFL)'); ylabel('Normalized Density');
title('TetR inverter plus/minus');
outputfig(h','TetR_induction_plusminus','plots');

TetR_SNR = SNR(results{1}.means(2),results{end}.means(2),results{1}.stds(2),results{end}.stds(2))
%    -0.0132 dB

%%%%%%%%%%%%%%%%%%%%
% Example of geometric vs. arithmetic statistics
arithmetic_mean = wmean(results{1}.bincenters,results{1}.bincounts(:,2)');
arithmetic_std = wstd(results{1}.bincenters,results{1}.bincounts(:,2)');
geometric_mean = geomean(results{1}.bincenters,results{1}.bincounts(:,2)');
geometric_std = geostd(results{1}.bincenters,results{1}.bincounts(:,2)');

% compute distribution functions
arithmetic_distribution = gmdistribution(arithmetic_mean,arithmetic_std^2,1);
% Note: gmdistribution's second argument is variance (i.e., std.dev.^2)
geometric_distribution = gmdistribution(log10(geometric_mean),log10(geometric_std)^2,1);
arithmetic_pdf = pdf(arithmetic_distribution,results{1}.bincenters');
geometric_pdf = pdf(geometric_distribution,log10(results{1}.bincenters'));

h = figure('PaperPosition',[1 1 6 4]);
semilogx(results{1}.bincenters,results{1}.bincounts(:,2)./max(results{1}.bincounts(:,2)),'k-'); hold on;
plot(results{1}.bincenters,geometric_pdf./max(geometric_pdf),'b--');
plot(results{1}.bincenters,arithmetic_pdf./max(arithmetic_pdf),'r--');
legend('TetR High','Geometric Fit','Arithmetic Fit');
xlim([1e0 1e4]); ylim([0 1.1]);
xlabel('RFP (MEFL)'); ylabel('Normalized Density');
title('TetR statistics comparison');
outputfig(h','statistics_comparison_TetR','plots');

