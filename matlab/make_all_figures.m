% MAKE_ALL_FIGURES  Reproduce the STORM synapse figures.
%
% Generates, for both ages:
%   * Synapse-density dot plots with mean +/- S.D. (Fig 4B / Fig S6A style)
%   * Volume & Intensity violin + box plots for VGluT2, Bassoon, Homer1
%     (Fig 4D / Fig S6C style)
%
% Outputs (.eps vector + .png) are written to ../output/.
% Run from the matlab/ directory, or add it to the MATLAB path.

thisDir = fileparts(mfilename('fullpath'));
outDir  = fullfile(thisDir, '..', 'output');
if ~isfolder(outDir), mkdir(outDir); end

data = load_storm_data();   % defaults to ../data/STORM_raw_spreadsheet.xlsx

% ----- Adult (Figure S6) -------------------------------------------------
plot_density(data, 'Adult', ...
    'SaveStem', fullfile(outDir, 'FigS6A_density_Adult'));
plot_volume_intensity(data, 'Adult', ...
    'SaveStem', fullfile(outDir, 'FigS6C_volint_Adult'));

% ----- P8 (Figure 4) -----------------------------------------------------
plot_density(data, 'P8', ...
    'SaveStem', fullfile(outDir, 'Fig4B_density_P8'));
plot_volume_intensity(data, 'P8', ...
    'SaveStem', fullfile(outDir, 'Fig4D_volint_P8'));

disp('All figures written to output/.');
