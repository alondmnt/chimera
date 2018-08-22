function [cars, cars_vec, err, filtered] = calc_cars(key, SA, ref, max_len, max_pos)
% [cars, cars_vec, err] = CALC_CARS(key, SA, ref, max_len, max_pos)
%  compute the ChimeraARS (Average Repeatetive Substring, Zur and Tuller
%  2014) for a given key.
%  an optimized implementation.
%
%   max_len: if provided, cARS will detect single substrings/blocks that
%       are larger than [max_len] and filter the entire gene of origin.
%   max_pos: if provided, cARS will detect genes that occur in a fraction
%       of positions that is larger than [max_pos] and filter them.
%
% Alon Diament / Tuller Lab, July 2015.

if nargin < 5 || max_pos <= 0 || ~isfinite(max_pos)
    max_pos = 1;
end
if nargin < 4 || max_len <= 0 || ~isfinite(max_len)
    max_len = length(key);
end

n = length(key);
cars_vec = zeros(n, 1);
cars_origin = zeros(n, 1);
err = false;
filtered = 0;

pos_list = true(n, 1);
while any(pos_list)
    pos = find(pos_list, 1);
    [substring, pid, homologs] = longest_prefix(key(pos:end), SA, ref, [], max_len);
    cars_origin(pos) = SA(pid, 2);
    cars_vec(pos) = length(substring);
    if isempty(substring)
        fprintf('empty substring at %d\n', pos);
        err = true;
    end
    if ~isempty(homologs)
        SA(homologs, 3) = 0;  % mask by 0 suffix frequency
        filtered = filtered + 1;
    end

    same = cars_origin == cars_origin(pos);
    if mean(same) > max_pos
        SA(SA(:, 2) == cars_origin(pos), 3) = 0;  % mask
        filtered = filtered + 1;
        pos_list(same) = true;
        cars_origin(same) = 0;
    else
        pos_list(pos) = false;
    end
end

cars = mean(cars_vec);
end
