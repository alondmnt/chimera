function [mapseq, B] = calc_cmap_posspec(key, SA, refAA, refNT, win_params)
% [mapseq, B] = CALC_CMAP_POSSPEC(key, SA, refAA, refNT, win_params)
%   compute the position-specific chimeraMap solution for a given key.
%   unlike the original chimeraMap (Zur and Tuller, 2014), blocks are
%   selected from windows in all reference sequences that are positioned
%   at the same distance (as in the target gene [key]) from ORF start/stop.
%
% Alon Diament / Tuller Lab, January 2018.

if ~isstruct(win_params)
    % assuming win_params is the size of the window
    win_params = struct('size', win_params, 'center', 0, ...
                        'by_start', true, 'by_stop', true, 'truncate_seq', false);
end

n = length(key);
B = cell(n, 3); % blocks
pos = 1; % position in key
mask = SA(:, 3);  % masking by suffix frequency
for blk = 1:n
    [SA, win_params] = select_window(SA, win_params, pos, pos-n-1);
    SA(:, 3) = min(mask, SA(:, 3));

    blockAA = longest_prefix(key(pos:end), SA, refAA, win_params);
    assert(~isempty(blockAA))
    [B{blk, 3}, B{blk, 1}, B{blk, 2}] = most_freq_prefix(blockAA, SA, refAA, refNT, win_params);
    pos = pos + length(blockAA);
    if pos > n
        break;
    end
end

B = B(1:blk, :);
mapseq = cat(2, B{:, 3});

assert(strcmp(nt2aa(mapseq, 'AlternativeStartCodons', false), key));
end


function [MF, gene, loc] = most_freq_prefix(pref, SA, refAA, refNT, win_params)
% finds the most frequent *NT* prefix in [SA] that codes the given 
% AA prefix [pref]. to this end, [SA] needs to be non-unique, so we can
% locate all occurences of the AA block in the genome.
COUNT_OVERELAPPING_WINS = false;  % simulating the parallel job version

nP = length(pref) - 1;
[~, left] = binary_search(pref, SA, refAA);
[~, right] = binary_search([pref, '~'], SA, refAA);
if left ~= right
    iSA = left : right - 1;
else
    % can happen at the end of SA
    iSA = left;
end
iSA(~SA(iSA, 3)) = [];  % 16/08/18: masking suffixes with 0 frequency
iSA(SA(iSA, 1) + nP > cellfun(@length, refAA(SA(iSA, 2)))) = [];  % 16/08/18: rudimentary check

% here we filter instead of truncating
count = zeros(length(iSA), 1);
if nargin > 4 && win_params.by_start
    if win_params.truncate_seq
        count = count + (SA(iSA, 1) >= win_params.win_start(1) & ...  % just in case (already handled in select_windows)
                         SA(iSA, 1) + nP <= win_params.win_start(2));  % entire seq within window bounds
    else
        count = count + (SA(iSA, 1) >= win_params.win_start(1) & ...
                         SA(iSA, 1) <= win_params.win_start(2));  % just in case (already handled in select_windows): seq starts within window bounds
    end
end
if nargin > 4 && win_params.by_stop
    if win_params.truncate_seq
        count = count + (SA(iSA, 4) >= win_params.win_stop(1) & ...
                         SA(iSA, 4) + nP <= win_params.win_stop(2));
    else
        count = count + (SA(iSA, 4) >= win_params.win_stop(1) & ...
                         SA(iSA, 4) <= win_params.win_stop(2));
    end
end
if any(count)
    iSA = iSA(count > 0);
    count = count(count > 0);
else
    assert(nargin < 5);  % we expect to have something if we got this far with windows given
end


if length(iSA) <= 1
    gene = SA(iSA, 2);  % left
    loc = SA(iSA, 1);  % left
    MF = refNT{gene}(3*(loc-1) + 1 : 3*(loc+nP));
    return;
end

all_blocks = arrayfun(@(g, pos) refNT{g}(3*(pos-1) + 1 : 3*(pos+nP)), ...
    SA(iSA, 2), SA(iSA, 1), 'UniformOutput', false);
nB = length(all_blocks);

[all_blocks, ind] = sort(all_blocks);
iSA = iSA(ind);
[all_blocks, iu] = unique(all_blocks, 'stable');
iSA = iSA(iu);
iu(end+1) = nB + 1;

if COUNT_OVERELAPPING_WINS && any(count)
    count = count(ind);
    freq = arrayfun(@(x, y) sum(count(x:y)), iu(1:end-1), iu(2:end)-1);
else
    freq = diff(iu);
end

[~, mostfreq] = max(freq);

MF = all_blocks{mostfreq}; % first in lexicographic order
gene = SA(iSA(mostfreq), 2);
loc = SA(iSA(mostfreq), 1);
end
