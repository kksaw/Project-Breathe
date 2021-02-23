function [xl, yl] = plotProbDistribution(ax, max_offset, pdoffset, xl, yl, marker, linewidth, markersize, markerec, markerfc)

% plotProbDistribution - plots a prob distribution
                        
line(ax, [0:max_offset-1], reshape(pdoffset, [max_offset,1]), ...
    'Color', markerec, ...
    'LineStyle', '-', ...
    'LineWidth', linewidth, ...
    'Marker', marker, ...
    'MarkerSize', markersize,...
    'MarkerEdgeColor', markerec,...
    'MarkerFaceColor', markerfc);
xl = [0 max(max_offset-1, xl(2))];
xlim(xl);
yl = [0 max(max(pdoffset) * 1.1, yl(2))];
ylim(yl);

end

