function [cars, cars_vec, err] = calc_cars(key, SA, ref)
% [cars, cars_vec, err] = CALC_CARS(key, SA, ref)
%  compute the ChimeraARS (Average Repeatetive Substring, Zur and Tuller
%  2014) for a given key.
%  an optimized implementation.
%
% Alon Diament, July 2015.

n = length(key);
cars_vec = zeros(1, n);
err = false;

for pos = 1:n
    substring = longest_prefix(key(pos:end), SA, ref);
    cars_vec(pos) = length(substring);
    if isempty(substring)
        fprintf('empty substring at %d\n', pos);
        err = true;
    end
end

cars = mean(cars_vec);
end
