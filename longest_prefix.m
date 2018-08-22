function [pref, pind, homologs] = longest_prefix(key, SA, ref, win_params, max_len)
%[pref, pind] = LONGEST_PREFIX(key, SA, ref, win_params, max_len)
%  find the longest prefix between [key] and the suffix array [SA],
%  and return its first index. original strings are stored in [ref] and 
%  their IDs referred to in [SA].
%
% Alon Diament, June 2015.
%
% changelog:
%   January 2018: [win_params] used to truncate suffixes (for PScMap-1).
%   August  2018: suffix masking.

where = binary_search(key, SA, ref);
where = min(max(1, where), size(SA, 1));

if nargin < 5
    max_len = length(key);
end
homologs = [];

n = Inf;  % prefix length
while n > max_len
    if isfinite(n)  % mask previous attempt (homolog suspect) and repeat
        if isempty(homologs)
            homologs = SA(:, 2) == SA(pind, 2);
        else
            homologs = homologs | (SA(:, 2) == SA(pind, 2));
        end
        SA(homologs, 3) = 0;
    end
    if nargin < 4 || isempty(win_params) || ~win_params.truncate_seq
        % when suffixes aren't truncated, the two adjacent suffixes contain the
        % longest common prefix
        nei = get_neighbors(where, SA);
        neighbors = arrayfun(@(x) {ref{SA(x, 2)}(SA(x, 1) : end)}, nei);
        [n, pind] = max(cellfun(@(x) count_common(key, x), neighbors));
        pind = nei(pind); % always return the first index that matched
        pref = key(1 : n);

    else
        [n, pind] = longest_trunc_prefix(key, SA, ref, where, win_params);
        pref = key(1 : n);
    end
end

end


function [n, pind] = longest_trunc_prefix(key, SA, ref, ind, win_params)
pind = ind;
n = 0;
iup = ind + 1;
idown = ind;

candidates_exist = true;
while candidates_exist
    % the main realization here is that we must continue to test neighbors
    % until the longest common prefix (LCP) between adjacent
    % (non-truncated) suffixes decreases below the longest truncated prefix
    % seen so far.
    % the stop condition checks the LCP of [iup] (new candidate) with
    % [pind] (current best) and the LCP of [pind] with [idown] (new
    % candidate).
    candidates_exist = false;

    iup = mask_suffix(SA, iup-1, -1, 0, Inf);
    if (iup > 0) && (count_common(ref{SA(iup, 2)}(SA(iup, 1):end), ...
                                  ref{SA(pind, 2)}(SA(pind, 1):end)) >= n)
        candidates_exist = true;
        tn = count_common(key, get_trunc_suffix(SA, ref, iup, win_params));
        if tn > n
            pind = iup;
            n = tn;
        end
    end

    idown = mask_suffix(SA, idown+1, +1, -Inf, size(SA, 1)+1);
    if (idown <= size(SA, 1)) && (count_common(ref{SA(idown, 2)}(SA(idown, 1):end), ...
                                               ref{SA(pind, 2)}(SA(pind, 1):end)) >= n)
        candidates_exist = true;
        tn = count_common(key, get_trunc_suffix(SA, ref, idown, win_params));
        if tn > n
            pind = idown;
            n = tn;
        end
    end
end
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


function neis = get_neighbors(loc, SA)
% it would seem that in order for masking to work properly we need to
% look both up (new) and down (legacy) the found suffix.

lo = mask_suffix(SA, max(1, loc-1), -1, 1, Inf);  % next lower neighbor
hi = mask_suffix(SA, min(size(SA, 1), loc+1), +1, -Inf, size(SA, 1));  % next upper neighbor

neis = [lo, loc, hi];
neis(~SA(neis, 3)) = [];
end
