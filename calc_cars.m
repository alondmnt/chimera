function [cars, cars_vec] = calc_cars(key, SA, ref)
% [cars, cars_vec] = CALC_CARS(key, SA, ref)
%  compute the ChimeraARS (Average Repeatetive Substring, Zur and Tuller
%  2014) for a given key.
%  an optimized implementation.
%
% Alon Diament, July 2015.

n = length(key);
cars_vec = zeros(1, n);
for pos = 1:n
    substring = longest_prefix(key(pos:end), SA, ref);
    cars_vec(pos) = length(substring);
end

cars = mean(cars_vec);
end
