function data = load_storm_data(xlsxPath)
% LOAD_STORM_DATA  Read the STORM raw spreadsheet into a struct of tables.
%
%   data = load_storm_data()            % uses ../data/STORM_raw_spreadsheet.xlsx
%   data = load_storm_data(xlsxPath)    % explicit path to the workbook
%
% Returns a struct with fields:
%   .density   table: Sample_ID, Age, Genotype, Replicate, Channel, Density
%   .VGluT2    table: Sample_ID, Age, Genotype, Replicate, Volume, Intensity
%   .Bassoon   same schema as .VGluT2
%   .Homer     same schema as .VGluT2
%
% Age is normalised to the categories "P8" and "Adult"; Genotype to
% "Het" and "KO". Volume and Intensity are the log-scaled values already
% present in the workbook (log(um^3) and log(A.U.) respectively).

    if nargin < 1 || isempty(xlsxPath)
        thisDir  = fileparts(mfilename('fullpath'));
        xlsxPath = fullfile(thisDir, '..', 'data', 'STORM_raw_spreadsheet.xlsx');
    end
    assert(isfile(xlsxPath), 'Workbook not found: %s', xlsxPath);

    % Preserve the original spreadsheet header text so column access does
    % not depend on MATLAB's version-specific name sanitisation. Columns are
    % then referenced by their exact header via getCol() below.

    % --- Synapse density -------------------------------------------------
    T = readSheet(xlsxPath, 'Synapse_density');
    dens = table;
    dens.Sample_ID = string(getCol(T, 'Sample_ID'));
    dens.Age       = normalizeAge(getCol(T, 'Age'));
    dens.Genotype  = normalizeGeno(getCol(T, 'Genotype'));
    dens.Replicate = double(getCol(T, 'Biological replicate'));
    dens.Channel   = string(getCol(T, 'Channel'));
    dens.Density   = double(getCol(T, 'Density'));
    dens = dens(~isnan(dens.Density), :);
    data.density = dens;

    % --- Volume / intensity sheets --------------------------------------
    sheets = struct('VGluT2','VGluT2_volume_intensity', ...
                    'Bassoon','Bassoon_volume_intensity', ...
                    'Homer','Homer_volume_intensity');
    prot = fieldnames(sheets);
    for i = 1:numel(prot)
        T = readSheet(xlsxPath, sheets.(prot{i}));
        vi = table;
        vi.Sample_ID = string(getCol(T, 'Sample_ID'));
        vi.Age       = normalizeAge(getCol(T, 'Age'));
        vi.Genotype  = normalizeGeno(getCol(T, 'Genotype'));
        vi.Replicate = double(getCol(T, 'Biological replicate'));
        vi.Volume    = double(getCol(T, 'Volume log(um3)'));
        vi.Intensity = double(getCol(T, 'Intensity log(A.U.)'));
        vi = vi(~isnan(vi.Volume) & ~isnan(vi.Intensity), :);
        data.(prot{i}) = vi;
    end
end

function T = readSheet(xlsxPath, sheet)
% Read a sheet keeping the original header strings as VariableNames.
    opts = detectImportOptions(xlsxPath, 'Sheet', sheet, ...
                               'VariableNamingRule', 'preserve');
    T = readtable(xlsxPath, opts);
end

function col = getCol(T, headerText)
% Fetch a column by its original header, tolerant of MATLAB having
% sanitised the name (e.g. older releases without 'preserve').
    vn = string(T.Properties.VariableNames);
    idx = find(vn == string(headerText), 1);
    if isempty(idx)
        % fall back to a sanitised match
        target = regexprep(string(headerText), '[^A-Za-z0-9]+', '_');
        san    = regexprep(vn, '[^A-Za-z0-9]+', '_');
        idx = find(san == target, 1);
    end
    assert(~isempty(idx), 'Column "%s" not found. Available: %s', ...
           headerText, strjoin(vn, ', '));
    col = T.(T.Properties.VariableNames{idx});
end

function a = normalizeAge(raw)
    s = string(raw);
    a = strings(size(s));
    a(ismember(lower(s), ["p8","p08"]))            = "P8";
    a(ismember(lower(s), ["adult","p60","p 60"]))  = "Adult";
    a(a=="") = s(a=="");   % keep anything unexpected as-is
    a = categorical(a, ["P8","Adult"]);
end

function g = normalizeGeno(raw)
    s = string(raw);
    g = strings(size(s));
    g(ismember(lower(s), ["het","opn4-het","opn4cre/+"])) = "Het";
    g(ismember(lower(s), ["ko","opn4-ko","opn4cre/cre"])) = "KO";
    g(g=="") = s(g=="");
    g = categorical(g, ["Het","KO"]);
end
