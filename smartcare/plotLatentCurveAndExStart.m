function [xl, yl] = plotLatentCurveAndExStart(ax, days, meancurve, xl, yl, nmean, ex_start, offset, colour)

% plotLatentCurveAndExStart - plots the latent curve along with the
% predicted exacerbation start

% plot mean curve (actual in dotted line, smoothed in solid line)
line(ax, days, meancurve + nmean, ...
    'Color', colour, ...
    'LineStyle', ':', ...
    'LineWidth', 1);
line(ax, days, smooth(meancurve + nmean, 5), ...
    'Color', colour, ...
    'LineStyle', '-', ...
    'LineWidth', 1);
xl = [min(min(days), xl(1)) max(max(days), xl(2))];
xlim(xl);
yl = [min(min((meancurve + nmean) * 0.99), yl(1)) max(max((meancurve + nmean) * 1.01), yl(2))];
ylim(yl);

% plot vertical line for predicted exacerbation start
line(ax, [ex_start + offset ex_start + offset] , yl, ...
    'Color', colour, ...
    'LineStyle', '-', ...
    'LineWidth', 0.5);
xl = [min(min(ex_start + offset), xl(1)) max(max(ex_start + offset), xl(2))];
xlim(xl);

% plot short vertical line for average exacerbation start indicator
line(ax, [ex_start ex_start], [yl(1), yl(1) + ((yl(2)-yl(1)) * 0.1)], ...
    'Color', colour, ...
    'LineStyle', ':', ...
    'LineWidth', 0.5);
xl = [min(min(ex_start), xl(1)) max(max(ex_start), xl(2))];
xlim(xl);

end

