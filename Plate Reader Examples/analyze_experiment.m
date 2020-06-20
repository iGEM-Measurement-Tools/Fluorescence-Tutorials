%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin by loading and verifying the data

datafile = 'iGEM 2019 Plate Reader Fluorescence Calibration - Example.xlsx';
template = iGEM_2019_plate_reader_fluorescence();
results = ExcelTemplateExtraction.extract(datafile,template);
validate_plate_Abs600(results);
validate_plate_fluorescence(results);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The results includes arrays of converted values for the data

estimated_cells = results.experiment_particles;
MEFL_per_cell = results.experiment_MEx_per_particle;

% In this case, each column is a replicate set
constructs = {'Negative','Positive','J23101','J23106','J23117','J23100','J23104','J23116'};
n_constructs = numel(constructs);

media_replicates = estimated_cells(:,9);
negative_control_replicates = MEFL_per_cell(:,1);

% Convert fluorescence to net values by subtracting autofluorescence
autofluorescence_mean = geomean(negative_control_replicates);
net_MEFL_per_cell = MEFL_per_cell - autofluorescence_mean;

% Values too close to autofluorescence / media cannot be distinguished from background
autofluorescence_std = geostd(negative_control_replicates); % geometric, because cells dominate
indistinguishable_MEFL = autofluorescence_mean*(autofluorescence_std^2 - 1);
media_std = std(media_replicates); % arithmetic, because instrument error dominates
indistinguishable_cells = 2*media_std;


% Check that absolute negative and positive control values are reasonable
% and that both of these sets of cells have grown well.
NEGATIVE_MAX = 1e3;
POSITIVE_MIN = 1e4;
assert(autofluorescence_mean < NEGATIVE_MAX);
assert(geomean(net_MEFL_per_cell(:,2)) > POSITIVE_MIN);
assert(geomean(estimated_cells(:,1)) > 1e7);
assert(geomean(estimated_cells(:,2)) > 1e7);

% Plot the control fluorescence values vs. thresholds
mean_MEFL = log10(geomean([negative_control_replicates net_MEFL_per_cell(:,2)]));
std_MEFL = log10(geostd([negative_control_replicates net_MEFL_per_cell(:,2)]));
h = figure('PaperPosition',[1 1 3 4]);
barwitherr([std_MEFL; std_MEFL]', mean_MEFL'); hold on;
plot([0 3],log10([NEGATIVE_MAX NEGATIVE_MAX]),'r--');
plot([0 3],log10([POSITIVE_MIN POSITIVE_MIN]),'r--');
xlabel('Construct');
xlim([0.5 2.5]); ylim([2 6]);
set(gca,'XTick',1:2);
set(gca,'XTickLabels',{'Raw Negative','Net Positive'});
ylabel('log_{10} Fluorescence (MEFL/cell)');
title('Control Values');
outputfig(h,'control_fluorescence','plots');

% Finally, plot the experimental results as bar graphs:
mean_MEFL = log10(geomean(net_MEFL_per_cell(:,3:n_constructs)));
std_MEFL = log10(geostd(net_MEFL_per_cell(:,3:n_constructs)));
h = figure('PaperPosition',[1 1 6 4]);
barwitherr([std_MEFL; std_MEFL]', mean_MEFL'); hold on;
plot([0 n_constructs+1],log10([indistinguishable_MEFL indistinguishable_MEFL]),'r--');
xlabel('Construct');
xlim([0.5 n_constructs-1.5]); ylim([2 6]);
set(gca,'XTick',1:(n_constructs-2));
set(gca,'XTickLabels',constructs(3:end));
ylabel('log_{10} Fluorescence (MEFL/cell)');
title('Construct Fluorescence');
outputfig(h,'construct_fluorescence','plots');

mean_cells = log10(geomean(estimated_cells(:,3:n_constructs)));
std_cells = log10(geostd(estimated_cells(:,3:n_constructs)));
h = figure('PaperPosition',[1 1 6 4]);
barwitherr([std_cells; std_cells]', mean_cells'); hold on;
plot([0 n_constructs+1],log10([indistinguishable_cells indistinguishable_cells]),'r--');
xlabel('Construct');
xlim([0.5 n_constructs-1.5]); ylim([6 9]);
set(gca,'XTick',1:(n_constructs-2));
set(gca,'XTickLabels',constructs(3:end));
ylabel('log_{10} Cell count');
title('Colony Size');
outputfig(h,'colony_size','plots');


