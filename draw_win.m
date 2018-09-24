function draw_win(h_axes, win_size, win_center, start_dist, stop_dist, truncate_seq)

if nargin == 0
    win_size = 70;
    win_center = 35;
    start_dist = true;
    stop_dist = true;
    truncate_seq = true;
    figure(1);
else
    axes(h_axes);
end
cla;

if ~win_size
    h_axes.XLim = [0, 100];
    h_axes.YLim = [0, 1];
    text(50, 0, 'not position-specific', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    return
end

if start_dist
    start_dist = 1.5*win_size;
end
if stop_dist
    stop_dist = 1.5*win_size;
end

pos = 0;
% L = start_dist + stop_dist;

win_x = pos - win_size/2 + win_center;
win_y = win_size / 4;
rectangle('Position', [win_x, 0, win_size, win_y], 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');
hold('on');

plot_arrow(win_x, 1.2*win_y, win_x + win_size, 1.2*win_y, 0.1);
plot_arrow(win_x + win_size, 1.2*win_y, win_x, 1.2*win_y, 0.1);
text(win_x + win_size/2, 1.2*win_y + win_size/16, 'window size', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

scatter(0, 0, 200, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
text(pos, -win_size/16, sprintf('position\nin target'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

if start_dist
    plot_arrow(-start_dist, 0, -0.07*start_dist, 0, 0.1);
    text(-start_dist/2, -win_size/16, sprintf('dist from\nSTART'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end
if stop_dist
    plot_arrow(stop_dist, 0, 0.07*stop_dist, 0, 0.1);
    text(stop_dist/2, -win_size/16, sprintf('dist from\nSTOP'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

axis('equal');
h_axes.XLim = [min(-1.5*win_size, win_x), max(1.5*win_size, win_x + win_size)];

str_block = repmat('ACGTACGTACGT', 1, 10);
for i = 1:100
    htext = text(win_x + win_size/3, win_y/2, str_block(1:i), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontName', 'fixedwidth', 'FontWeight', 'bold', 'FontSize', 12);
    if truncate_seq && (htext.Extent(1) + htext.Extent(3) > win_x + win_size)
        delete(htext);
        break
    end
    delete(htext);
end
text(win_x + win_size/3, win_y/2, str_block(1:i-1), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontName', 'fixedwidth', 'FontWeight', 'bold', 'FontSize', 12);

hold('off');

h_axes.XLim = [min(-1.5*win_size, win_x), max(1.5*win_size, win_x + win_size)];
set(gca, 'XTick', [], 'YTick', []);


function plot_arrow(x1, y1, x2, y2, head_size)
line([x1, x2], [y1, y2], 'Color', 'k', 'LineWidth', 2);
m = [x2 - x1, y2 - y1];

head_base = [x1, y1] + (1-head_size) * m;
head1 = head_base + 0.5*head_size * [0, m(1)];
head2 = head_base - 0.5*head_size * [0, m(1)];

line([head1(1), x2], [head1(2), y2], 'Color', 'k', 'LineWidth', 2);
line([head2(1), x2], [head2(2), y2], 'Color', 'k', 'LineWidth', 2);
