function [exists, where] = binary_search(key, SA, ref)
%[exists, where] = BINARY_SEARCH(key, S, ref)
%  find whether string [key] exists in sorted array [SA] and: 
%  (1) return its index; or 
%  (2) return before which index it should be inserted to preserve sorting.
%  original strings should be given in [ref] and referred to in [SA].
%
% currently a [key] that should be place *last* in the SA will return the 
% last index (instead of last + 1). we use that to our advantage.
%
% Alon Diament, July 2015.

nS = size(SA, 1);
low = 1;
high = nS;
i = 0;
while low < high
    i = i + 1;
    mid = floor((low + high) / 2);
    mid = mask_suffix(SA, mid, -1, low, Inf);

    if is_greater(key, mid, SA, ref)
        low = mask_suffix(SA, mid+1, +1, -Inf, high);
    else
        high = mask_suffix(SA, mid, -1, low, Inf);
    end
end
where = low;

if where == 1
    where = mask_suffix(SA, where, +1, -Inf, nS);
elseif where == nS
    where = mask_suffix(SA, where, -1, 1, Inf);
end
assert(SA(where, 3) > 0);

if where == nS  % case: end of SA
    if is_greater(key, where, SA, ref)
        where = nS + 1;
    end
end
if where > nS
    exists = false;
    return;
end

if strcmp(key, ref{SA(where, 2)}(SA(where, 1) : end))
    exists = true;
else
    exists = false; % insert [key] before [where]
end


function res = is_greater(key, ind, SA, ref)
[~, idx] = sort({key, ref{SA(ind, 2)}(SA(ind, 1) : end)});
res = idx(1) > idx(2);
