function [SA, win_params, empty_SA] = select_window(SA, win_params, pos_start, pos_stop, mask)
% [SA, win_params, empty_SA] = SELECT_WINDOW(SA, win_params, pos_start, pos_stop, mask)
% filter only relevant suffixes, i.e.:
% 1. starting at the same distance from start as [start_pos], within
% [max_dist] codons; or:
% 2. starting with at the same distance from stop as [end_pos], within
% [max_dist] codons.
%
% Alon Diament, Tuller Lab, Januray 2018.

do_binary = size(SA, 2) > 5;

max_dist = win_params.size / 2;
SA(:, 3) = 0;  % masking suffixes using 0 frequency
empty_SA = true;

if win_params.by_start
    win_start = [pos_start - max_dist + win_params.center, ...
                 pos_start + max_dist + win_params.center - 1];  % relative to start codon
    win_start = ceil(win_start);  % handling odd numbers

    if ~do_binary
        valid = (win_start(1) <= SA(:, 1) & SA(:, 1) <= win_start(2));
    else
        valid = binary_search(SA, 1, 5, win_start(1), win_start(2));
    end
    valid(mask(valid)) = [];
    SA(valid, 3) = 1;
    empty_SA = empty_SA & isempty(valid);
    win_params.win_start = win_start;
end

if win_params.by_stop
    win_stop(2) = pos_stop + max_dist + win_params.center -1;  % relative to stop codon
    win_stop(1) = pos_stop - max_dist + win_params.center;
    win_stop = ceil(win_stop);

    if ~do_binary
        valid = (win_stop(1) <= SA(:, 4) & SA(:, 4) <= win_stop(2));
    else
        valid = binary_search(SA, 4, 6, win_stop(1), win_stop(2));
    end
    valid(mask(valid)) = [];
    SA(valid, 3) = 1;
    empty_SA = empty_SA & isempty(valid);
    win_params.win_stop = win_stop;
end

end


function valid = binary_search(SA, val_col, idx_col, start, stop)
assert(start <= stop);

low = 1;
high = size(SA, 1);
while low < high
    mid = floor((low + high) / 2);

    if SA(SA(mid, idx_col), val_col) < start
        low = mid + 1;
    else
        high = mid;
    end
end
i = low;

high = size(SA, 1);
while low < high
    mid = floor((low + high) / 2);
    if SA(SA(mid, idx_col), val_col) > stop
        high = mid;
    else
        low = mid + 1;
    end
end
if stop == SA(SA(end, idx_col), val_col)
    j = low;
else
    j = low - 1;
end

valid = SA(i:j, idx_col);
end
