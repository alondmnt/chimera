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

low = 1;
high = size(SA, 1);
while low < high
    mid = floor((low + high) / 2);
    
    [~, idx] = sort({key, ref{SA(mid, 2)}(SA(mid, 1) : end)});
    isGreater = idx(1) > idx(2);
    if isGreater
        low = mid + 1;
    else
        high = mid;
    end
end
where = low;

if where == size(SA, 1)  % case: end of SA
    [~, idx] = sort({key, ref{SA(end, 2)}(SA(end, 1) : end)});
    isGreater = idx(1) > idx(2);
    if isGreater
        where = size(SA, 1) + 1;
        exists = false;
        return
    end
end

if where > size(SA, 1)
    exists = false;
    return;
end
% low == high
if strcmp(key, ref{SA(low, 2)}(SA(low, 1) : end))
    exists = true;
else
    exists = false; % insert [key] before [where]
end
