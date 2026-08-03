function fig = plot_volume_intensity(data, age, varargin)
% PLOT_VOLUME_INTENSITY  Violin + box-and-whisker panels for the three
% synaptic proteins (VGluT2, Bassoon, Homer1), showing Volume (left column)
% and Intensity (right column) for Het (black) vs KO (red).
%
%   fig = plot_volume_intensity(data, 'Adult')   % Figure S6C style
%   fig = plot_volume_intensity(data, 'P8')      % Figure 4D style
%
% Each violin is a kernel-density outline; the overlaid box shows the
% median and interquartile range, with whiskers extending to 1.5x IQR.
%
% Name/value options:
%   'SaveStem'  path stem; if given, saves <stem>.eps and <stem>.png

    p = inputParser;
    addParameter(p, 'SaveStem', "");
    parse(p, varargin{:});
    saveStem = string(p.Results.SaveStem);

    age      = string(age);
    proteins = ["VGluT2","Bassoon","Homer"];
    protLbl  = ["VGluT2","Bassoon","Homer1"];
    measures = ["Volume","Intensity"];
    measLbl  = ["Volume log(\mum^3)","Intensity log(A.U.)"];
    % Beeswarm point-cloud colours: Het/WT gray, KO red.
    cloudColors = struct('Het', [0.55 0.55 0.55], 'KO', [0.85 0.10 0.10]);
    genos    = ["Het","KO"];

    fig = figure('Color','w','Units','inches','Position',[1 1 9 3.4]);
    tl  = tiledlayout(fig, 1, 6, 'Padding','compact', 'TileSpacing','compact');

    for pi = 1:numel(proteins)
        T = data.(char(proteins(pi)));
        T = T(T.Age == age, :);
        for mi = 1:numel(measures)
            ax = nexttile(tl); hold(ax,'on');
            allv = T.(char(measures(mi)));
            panelSpan = max(allv) - min(allv);
            for k = 1:numel(genos)
                g    = genos(k);
                mask = T.Genotype == g;
                v    = T.(char(measures(mi)))(mask);
                reps = T.Replicate(mask);
                % per-animal (per-replicate) means, N = 6
                ur   = unique(reps);
                am   = arrayfun(@(r) mean(v(reps==r)), ur);
                drawViolinBox(ax, k, v, cloudColors.(char(g)), am, panelSpan);
            end
            xlim(ax, [0.4 2.6]);
            set(ax, 'XTick', [1 2], 'XTickLabel', {'Het','KO'});
            ylabel(ax, measLbl(mi));
            title(ax, protLbl(pi), 'FontWeight','normal');
            ax.FontName = 'Arial'; ax.FontSize = 8;
            ax.Box = 'off'; ax.TickDir = 'out';
        end
    end
    title(tl, sprintf('%s synaptic protein distributions', age), ...
          'FontName','Arial','FontWeight','normal');

    if strlength(saveStem) > 0
        exportgraphics(fig, saveStem + ".png", 'Resolution', 300);
        % PDF preserves the point-cloud transparency (EPS does not).
        exportgraphics(fig, saveStem + ".pdf", 'ContentType', 'vector');
    end
end

function drawViolinBox(ax, xc, v, cloudC, animalMeans, yspan)
% Every synapse is drawn as a circle; the horizontal spread of each point is
% scaled by the local point density so the cloud itself forms the violin
% silhouette. The point cloud is drawn in cloudC (gray for Het, red for KO);
% the per-animal means and the box-and-whisker are always black. The N = 6
% mean dots are fanned out horizontally where they would otherwise overlap.
    v = v(~isnan(v));
    if numel(v) < 2, return; end
    width = 0.36;
    black = [0 0 0];

    % --- local density -> per-point horizontal offset (violin from points) ---
    if all(v == v(1))
        off = zeros(size(v));
    else
        dens = ksdensity(v, v);          % density evaluated at each point
        dens = dens / max(dens);         % 0..1
        rng(0);                          % reproducible cloud
        off = (rand(numel(v),1)*2 - 1) .* dens(:) * width;
    end
    scatter(ax, xc + off(:), v(:), 1.5, cloudC, 'filled', ...
            'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');

    % --- per-animal means (N = 6) overlaid as black circles, de-overlapped ---
    if nargin >= 5 && ~isempty(animalMeans)
        am = animalMeans(~isnan(animalMeans));
        if nargin < 6 || isempty(yspan) || yspan == 0
            yspan = max(am) - min(am) + eps;
        end
        dx = dodgeOverlaps(am, 0.035*yspan, 0.05);
        scatter(ax, xc + dx(:), am(:), 16, black, 'filled', ...
                'MarkerEdgeColor', 'w', 'LineWidth', 0.4);
    end

    % --- box-and-whisker statistics (black lines) ---
    q1 = quantile(v, 0.25);
    q2 = median(v);
    q3 = quantile(v, 0.75);
    iqr = q3 - q1;
    loW = max(min(v), q1 - 1.5*iqr);
    hiW = min(max(v), q3 + 1.5*iqr);

    bw = 0.12;   % box half-width
    % whiskers
    plot(ax, [xc xc], [loW q1], 'Color', black, 'LineWidth', 1.0);
    plot(ax, [xc xc], [q3 hiW], 'Color', black, 'LineWidth', 1.0);
    plot(ax, xc+[-bw bw], [loW loW], 'Color', black, 'LineWidth', 1.0);
    plot(ax, xc+[-bw bw], [hiW hiW], 'Color', black, 'LineWidth', 1.0);
    % box (unfilled so the point cloud shows through)
    rectangle(ax, 'Position', [xc-bw, q1, 2*bw, q3-q1], ...
              'EdgeColor', black, 'FaceColor', 'none', 'LineWidth', 1.2);
    % median
    plot(ax, xc+[-bw bw], [q2 q2], 'Color', black, 'LineWidth', 1.6);
end
