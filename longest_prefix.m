function [pref, pind] = longest_prefix(key, SA, ref, win_params)
%[pref, pind] = LONGEST_PREFIX(key, SA, ref, max_pos_start, max_pos_end)
%  find the longest prefix between [key] and the suffix array [SA],
%  and return its first index. original strings are stored in [ref] and 
%  their IDs referred to in [SA].
%
% Alon Diament, June 2015.
%
% changelog:
%   January 2018: [win_params] used to truncate suffixes (for PScMap-1).

[exists, where] = binary_search(key, SA, ref);
if exists
    pind = where;
    if nargin < 4
        pref = key;
    else
        pref = get_trunc_suffix(SA, ref, where, win_params);
    end
    return;
end

where = min(max(2, where), size(SA, 1));
if nargin < 4 || ~win_params.truncate_seq
    % when suffixes aren't truncated, the two adjacent suffixes contain the
    % longest common prefix
    neighbors = {ref{SA(where-1, 2)}(SA(where-1, 1) : end); ...
                 ref{SA(where, 2)}(SA(where, 1) : end)};
    [n, pind] = max(cellfun(@(x) count_common(key, x), neighbors));
    pind = where + pind - 2; % always return the first index that matched
    pref = key(1 : n);
else
    [n, pind] = longest_trunc_prefix(key, SA, ref, where, win_params);
    pref = key(1 : n);
end

end


function [n, pind] = longest_trunc_prefix(key, SA, ref, ind, win_params)
i = 0;
pind = ind;
n = 0;

candidates_exist = true;
while candidates_exist
    % (SA(ind-1-i, LCP_IND) > n) || (SA(ind-1+i, LCP_IND) > n)
    % the main realization here is that we must continue to test neighbors
    % until the longest common prefix (LCP) between adjacent
    % (non-truncated) suffixes decreases below the longest truncated prefix
    % seen so far.
    % we check the LCP of [ind-1-i] (new candidate) with [ind-i] (prev cand)
    % and the LCP of [ind-1+i] (prev cand) with [ind+i] (new cand)
    candidates_exist = false;

    iup = ind - 1 - i;
    if (iup > 0) && (count_common(ref{SA(iup, 2)}(SA(iup, 1):end), ...
                                  ref{SA(pind, 2)}(SA(pind, 1):end)) >= n)
        candidates_exist = true;
        tn = count_common(key, get_trunc_suffix(SA, ref, iup, win_params));
        if tn > n
            pind = iup;
            n = tn;
        end
    end

    idown = ind + i;
    if (idown <= size(SA, 1)) && (count_common(ref{SA(idown, 2)}(SA(idown, 1):end), ...
                                               ref{SA(pind, 2)}(SA(pind, 1):end)) >= n)
        candidates_exist = true;
        tn = count_common(key, get_trunc_suffix(SA, ref, idown, win_params));
        if tn > n
            pind = idown;
            n = tn;
        end
    end
    i = i + 1;
end
% fprintf('%d steps (selected %d)\n', i, pind - ind);

end


function suf = get_trunc_suffix(SA, ref, ind, win_params)
% in order to reproduce the behavior of two legitimate windows (one
% measured from START and the other measured from STOP), need to take the
% maximal suffix that is included within any of the windows.

if nargin < 4
    n = length(ref{SA(ind, 2)}) - SA(ind, 1);  % all suffix
else
    n = 0;  % according to given windows
end

if nargin > 3 && win_params.by_start && (SA(ind, 1) >= win_params.win_start(1))
%     assert(all(win_params.win_start > 0), 'start position must be positive (relative to START)');
    nstart = win_params.win_start(2) - SA(ind, 1);
    if nstart >= 0
        n = max(n, nstart);
    end
end
if nargin > 3 && win_params.by_stop && (SA(ind, 4) >= win_params.win_stop(1))
%     assert(all(win_params.win_stop < 0), 'end position must be negative (relative to STOP)');
    nstop = win_params.win_stop(2) - SA(ind, 4);
    if nstop >= 0
        n = max(n, nstop);
    end
end

suf = ref{SA(ind, 2)}(SA(ind, 1) : min(end, SA(ind, 1) + n));
end


function c = count_common(s1, s2)
% counting the number of first matching characters
n = min(length(s1), length(s2));
c = find(s1(1:n) - s2(1:n), 1, 'first') - 1;
if isempty(c)
    c = n;
end
end
