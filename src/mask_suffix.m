function ind = mask_suffix(SA, ind, step, min_ind, max_ind)
% getting the next suffix (up/down) in SA for which the frequency is
% non-zero.

while (min_ind < ind && ind < max_ind) && ~SA(ind, 3)
    ind = ind + step;
end
