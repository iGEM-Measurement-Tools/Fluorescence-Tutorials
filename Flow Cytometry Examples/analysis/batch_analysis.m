% This template shows how to perform a simple batch analysis of a set of conditions
% Each color is analyzed independently
TASBEConfig.checkpoint(TASBEConfig.checkpoints());

% load the color model
load('../colormodel/CM140521.mat');
% can also add filters, such as gating out all "low transfection" red less than 10^6 MEFL:
%CM = add_postfilter(CM,RangeFilter('PE-Tx-Red-YG-A',[1e6 inf]));


% set up metadata
experimentName = 'TetR Transfer Curve';

% Configure the analysis
% Analyze on a histogram of 10^[first] to 10^[third] ERF, with bins every 10^[second]
bins = BinSequence(0,0.1,6,'log_bins');

% Designate which channels have which roles
AP = AnalysisParameters(bins,{});
% Ignore any bins with less than valid count as noise
AP=setMinValidCount(AP,100');
% Ignore any raw fluorescence values less than this threshold as too contaminated by instrument noise
AP=setPemDropThreshold(AP,5');
% Add autofluorescence back in after removing for compensation?
AP=setUseAutoFluorescence(AP,false');
% By default, analysis tries to fit constitutive to transformed and non-transformed components
% If your distribution is more complex or less complex, you can change the number of components
% AP=setNumGaussianComponents(AP,3);

% Make a map of condition names to file sets
stem0521 = '../20140521_Ara-TetR_Inv_fcs/Specimen_';
file_pairs = {...
  'Ara 1000 uM',    {DataFile('fcs', [stem0521 '001_A1_A01_049.fcs']), DataFile('fcs', [stem0521 '002_A2_A02_057.fcs']), DataFile('fcs', [stem0521 '003_A3_A03_065.fcs'])};
  'Ara  100 uM',    {DataFile('fcs', [stem0521 '001_B1_B01_050.fcs']), DataFile('fcs', [stem0521 '002_B2_B02_058.fcs']), DataFile('fcs', [stem0521 '003_B3_B03_066.fcs'])};
  'Ara   50 uM',    {DataFile('fcs', [stem0521 '001_C1_C01_051.fcs']), DataFile('fcs', [stem0521 '002_C2_C02_059.fcs']), DataFile('fcs', [stem0521 '003_C3_C03_067.fcs'])};
  'Ara   10 uM',    {DataFile('fcs', [stem0521 '001_D1_D01_052.fcs']), DataFile('fcs', [stem0521 '002_D2_D02_060.fcs']), DataFile('fcs', [stem0521 '003_D3_D03_068.fcs'])};
  'Ara    1 uM',    {DataFile('fcs', [stem0521 '001_E1_E01_053.fcs']), DataFile('fcs', [stem0521 '002_E2_E02_061.fcs']), DataFile('fcs', [stem0521 '003_E3_E03_069.fcs'])};
  'Ara    0.5 uM',  {DataFile('fcs', [stem0521 '001_F1_F01_054.fcs']), DataFile('fcs', [stem0521 '002_F2_F02_062.fcs']), DataFile('fcs', [stem0521 '003_F3_F03_070.fcs'])};
  'Ara    0.1 uM',  {DataFile('fcs', [stem0521 '001_G1_G01_055.fcs']), DataFile('fcs', [stem0521 '002_G2_G02_063.fcs']), DataFile('fcs', [stem0521 '003_G3_G03_071.fcs'])};
  'Ara    0 uM',    {DataFile('fcs', [stem0521 '001_H1_H01_056.fcs']), DataFile('fcs', [stem0521 '002_H2_H02_064.fcs']), DataFile('fcs', [stem0521 '003_H3_H03_072.fcs'])};
  };

n_conditions = size(file_pairs,1);

% Create point cloud files
% TASBEConfig.set('flow.outputPointCloud','true');
% TASBEConfig.set('flow.pointCloudPath','csv/');
% TASBEConfig.set('flow.pointCloudFileType', 1);

% Execute the actual analysis
TASBEConfig.set('OutputSettings.StemName','Induction');
TASBEConfig.set('OutputSettings.FixedBinningAxis',[1e0 1e6]);
[results, sampleresults] = per_color_constitutive_analysis(CM,file_pairs,{'GFP','RFP'},AP);

% Make output plots
plot_batch_histograms(results,sampleresults,CM); % linespecs obtained from CM
% can enter own linespecs for plot_batch_histograms:
% plot_batch_histograms(results,sampleresults,CM,{'b','g','r'});

[statisticsFile, histogramFile] = serializeBatchOutput(file_pairs, CM, AP, sampleresults);

save('TetR-batch.mat','AP','bins','file_pairs','results','sampleresults');
