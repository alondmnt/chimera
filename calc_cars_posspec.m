function [cars, cars_vec, err, filtered] = calc_cars_posspec(key, SA, ref, win_params, max_len, max_pos)
% [cars, cars_vec, err, filtered] = CALC_CARS_POSSPEC(key, SA, ref, win_params, max_len, max_pos)
%   compute the position-specific chimeraARS score for a given key.
%   unlike the original chimeraARS (Zur and Tuller, 2014), blocks are
%   selected from windows in all reference sequences that are positioned
%   at the same distance (as in the target gene [key]) from ORF start/stop.
%
%   max_len: if provided, cARS will detect single substrings/blocks that
%       are larger than [max_len] and filter the entire gene of origin.
%   max_pos: if provided, cARS will detect genes that occur in a fraction
%       of positions that is larger than [max_pos] and filter them.
%
% Alon Diament / Tuller Lab, June 2018.

if ~isstruct(win_params)
    % assuming win_params is the size of the window
    win_params = struct('size', win_params, 'center', 0, ...
                        'by_start', true, 'by_stop', true, 'truncate_seq', false);
end
if nargin < 6 || max_pos <= 0 || ~isfinite(max_pos)
    max_pos = 1;
end
if nargin < 5 || max_len <= 0 || ~isfinite(max_len)
    max_len = length(key);
end

n = length(key);
cars_vec = zeros(n, 1);
cars_origin = zeros(n, 1);
mask = ~SA(:, 3);  % masking by suffix frequency
err = false(2, 1);
filtered = 0;

pos_list = true(n, 1);
while any(pos_list)
    pos = find(pos_list, 1);
    [SA, win_params, empty_SA] = select_window(SA, win_params, pos, pos-n-1, mask);
    if empty_SA
        fprintf('empty window at %d\n', pos);
        cars_vec(pos) = 0;
        cars_origin(pos) = NaN;
        pos_list(pos) = false;
        err(2) = true;
        continue
    end

    [substring, pid, homologs] = longest_prefix(key(pos:end), SA, ref, win_params, max_len);
    cars_origin(pos) = SA(pid, 2);
    cars_vec(pos) = length(substring);
    if isempty(substring)
        fprintf('empty substring at %d\n', pos);
        err(1) = true;
    end
    if ~isempty(homologs)
        mask(homologs) = true;
        filtered = filtered + 1;
    end

    same = cars_origin == cars_origin(pos);
    if mean(same) > max_pos
        mask(SA(:, 2) == cars_origin(pos)) = true;
        filtered = filtered + 1;
        pos_list(same) = true;
        cars_origin(same) = 0;
    else
        pos_list(pos) = false;
    end
end

cars = mean(cars_vec);
end
