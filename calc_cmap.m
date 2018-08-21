function [mapseq, B, err] = calc_cmap(key, SA, refAA, refNT)
% [mapseq, B, err] = CALC_CMAP(key, SA, refAA, refNT)
%  compute the ChimeraMap (Zur and Tuller 2014) solution for a given key.
%  an optimized implementation.
%
% Alon Diament, July 2015.

n = length(key);
B = cell(n, 3); % blocks
pos = 1; % position in key
err = false(3, 1);

for blk = 1:n
    blockAA = longest_prefix(key(pos:end), SA, refAA);
    if isempty(blockAA)
        fprintf('empty block at %d\n', pos);
        B(blk, :) = {NaN, NaN, '---'};
        pos = pos + 1;
        err(1) = true;
        continue
    end

    [B{blk, 3}, B{blk, 1}, B{blk, 2}] = most_freq_prefix(blockAA, SA, refAA, refNT);
    pos = pos + length(blockAA);
    if pos > n
        break
    end
end

B = B(1:blk, :);
mapseq = cat(2, B{:, 3});

if ~any(err)
	err(3) = ~strcmp(nt2aa(mapseq, 'AlternativeStartCodons', false), key);
end
end


function [MF, gene, loc] = most_freq_prefix(pref, SA, refAA, refNT)
% finds the most frequent *NT* prefix in [SA] that codes the given 
% AA prefix [pref]. to this end, [SA] needs to be non-unique, so we can
% locate all occurences of the AA block in the genome.

nP = length(pref) - 1;
left = binary_search(pref, SA, refAA);
right = binary_search([pref, '~'], SA, refAA);
iSA = left : right - 1;
iSA(~SA(iSA, 3)) = [];  % 16/08/18: masking suffixes with 0 frequency
iSA(SA(iSA, 1) + nP > cellfun(@length, refAA(SA(iSA, 2)))) = [];  % 16/08/18: rudimentary check

if length(iSA) <= 1
    gene = SA(left, 2);
    loc = SA(left, 1);
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
freq = diff(iu);
[~, mostfreq] = max(freq);

MF = all_blocks{mostfreq}; % first in lexicographic order
gene = SA(iSA(mostfreq), 2);
loc = SA(iSA(mostfreq), 1);

end
