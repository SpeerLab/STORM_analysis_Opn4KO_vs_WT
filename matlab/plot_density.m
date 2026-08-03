function fig = plot_density(data, age, varargin)
% PLOT_DENSITY  Synapse-density dot plot with mean +/- S.D. error bars.
%
%   fig = plot_density(data, 'Adult')   % Figure S6A style
%   fig = plot_density(data, 'P8')      % Figure 4B style
%
% Each biological replicate (N = 6 per genotype) is plotted as an
% individual jittered dot: Het in black, KO in red. A horizontal bar marks
% the group mean and the error bars span +/- 1 standard deviation.
%
% Name/value options:
%   'Channel'   channel to plot from the density sheet (default "VGluT2")
%   'SaveStem'  path stem; if given, saves <stem>.eps and <stem>.png

    p = inputParser;
    addParameter(p, 'Channel', "VGluT2");
    addParameter(p, 'SaveStem', "");
    parse(p, varargin{:});
    channel  = string(p.Results.Channel);
    saveStem = string(p.Results.SaveStem);

    age = string(age);
    colors = struct('Het', [0 0 0], 'KO', [0.85 0.10 0.10]);
    genos  = ["Het","KO"];

    D = data.density;
    D = D(D.Age == age & D.Channel == channel, :);

    fig = figure('Color','w','Units','inches','Position',[1 1 2.2 3.0]);
    ax = axes(fig); hold(ax,'on');

    xpos = [1 2];
    yrange = max(D.Density) * 1.15;   % matches ylim below
    for k = 1:numel(genos)
        g  = genos(k);
        v  = D.Density(D.Genotype == g);
        c  = colors.(char(g));
        mu = mean(v);
        sd = std(v, 0);   % sample S.D. (N-1)

        % individual animals near the error bar; nudge apart if they overlap
        dx = dodgeOverlaps(v, 0.03*yrange, 0.055);
        scatter(ax, xpos(k)+dx(:), v, 30, 'MarkerFaceColor', c, ...
                'MarkerEdgeColor', c, 'MarkerFaceAlpha', 0.9);

        % mean bar + S.D. error bars
        errorbar(ax, xpos(k), mu, sd, 'Color', c, 'LineWidth', 1.2, ...
                 'CapSize', 12, 'LineStyle','none');
        plot(ax, xpos(k)+[-0.22 0.22], [mu mu], 'Color', c, 'LineWidth', 1.6);
    end

    xlim(ax, [0.4 2.6]);
    set(ax, 'XTick', xpos, 'XTickLabel', {'Het','KO'});
    ylabel(ax, 'Synapse density (per \mum^3)');
    title(ax, sprintf('%s  (%s)', channel, age), 'FontWeight','normal');
    ax.FontName = 'Arial'; ax.FontSize = 9;
    ax.Box = 'off'; ax.TickDir = 'out';
    ylim(ax, [0 max(D.Density)*1.15]);

    if strlength(saveStem) > 0
        exportgraphics(fig, saveStem + ".png", 'Resolution', 300);
        % PDF preserves the marker transparency (EPS does not).
        exportgraphics(fig, saveStem + ".pdf", 'ContentType', 'vector');
    end
end
