function [where, exists] = binary_search(key, SA, ref)
%[where, exists] = BINARY_SEARCH(key, S, ref)
%  find whether string [key] exists in sorted array [SA] and: 
%  (1) return its index; or 
%  (2) return before which index it should be inserted to preserve sorting.
%  original strings should be given in [ref] and referred to in [SA].
%
% Alon Diament, July 2015.

nS = size(SA, 1);
low = 1;
high = nS;
i = 0;
while low < high
    i = i + 1;
    mid = floor((low + high) / 2);

    if is_greater(key, mid, SA, ref)
        low = mid + 1;
    else
        high = mid;
    end
end
where = low;

if where < nS
    where = mask_suffix(SA, where, +1, -Inf, nS);
end
if where == nS  % case: end of SA
    where = mask_suffix(SA, where, -1, 1, Inf);  % last unmasked suffix
    if is_greater(key, where, SA, ref)
        where = where + 1;
    end
end

if nargout == 1
    return
end

if where > nS
    exists = false;
elseif strcmp(key, ref{SA(where, 2)}(SA(where, 1) : end))
    exists = true;
else
    exists = false; % insert [key] before [where]
end


function res = is_greater(key, ind, SA, ref)
[~, idx] = sort({key, ref{SA(ind, 2)}(SA(ind, 1) : end)});
res = idx(1) > idx(2);
