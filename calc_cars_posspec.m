function [cars, cars_vec, err] = calc_cars_posspec(key, SA, ref, win_params)
% [cars, cars_vec, err] = CALC_CARS_POSSPEC(key, SA, ref, win_params)
%   compute the position-specific chimeraARS score for a given key.
%   unlike the original chimeraARS (Zur and Tuller, 2014), blocks are
%   selected from windows in all reference sequences that are positioned
%   at the same distance (as in the target gene [key]) from ORF start/stop.
%
% Alon Diament / Tuller Lab, June 2018.

if ~isstruct(win_params)
    % assuming win_params is the size of the window
    win_params = struct('size', win_params, 'center', 0, ...
                        'by_start', true, 'by_stop', true, 'truncate_seq', false);
end

n = length(key);
cars_vec = zeros(1, n);
mask = ~SA(:, 3);  % masking by suffix frequency
err = false(2, 1);

for pos = 1:n
    [SA, win_params, empty_SA] = select_window(SA, win_params, pos, pos-n-1, mask);
    if empty_SA
        fprintf('empty window at %d\n', pos);
        cars_vec(pos) = 0;
        err(2) = true;
        continue
    end

    substring = longest_prefix(key(pos:end), SA, ref, win_params);
    cars_vec(pos) = length(substring);
    if isempty(substring)
        fprintf('empty substring at %d\n', pos);
        err(1) = true;
    end
end

cars = mean(cars_vec);
end
