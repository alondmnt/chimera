function [mapseq, B, err, filtered] = calc_cmap(key, SA, refAA, refNT, max_len, max_pos)
% [mapseq, B, err, filtered] = CALC_CMAP(key, SA, refAA, refNT, max_len, max_pos)
%  compute the ChimeraMap (Zur and Tuller 2014) solution for a given key.
%  an optimized implementation.
%
%   max_len: if provided, cARS will detect single substrings/blocks that
%       are larger than [max_len] and filter the entire gene of origin.
%   max_pos: if provided, cARS will detect genes that occur in a fraction
%       of positions that is larger than [max_pos] and filter them.
%
% Alon Diament / Tuller Lab, July 2015.

if nargin < 6 || max_pos <= 0 || ~isfinite(max_pos)
    max_pos = 1;
end
if nargin < 5 || max_len <= 0 || ~isfinite(max_len)
    max_len = length(key);
end

n = length(key);
B = cell(n, 3); % blocks
cmap_origin = zeros(n, 2);
err = false(3, 1);
filtered = 0;

pos = 1; % position in key
blk = 0;
while pos <= n
    blk = blk + 1;
    [blockAA, ~, homologs] = longest_prefix(key(pos:end), SA, refAA, [], max_len);
    m = length(blockAA);
    if isempty(blockAA)
        fprintf('empty block at %d\n', pos);
        B(blk, :) = {NaN, NaN, '---'};
        cmap_origin(pos, :) = [NaN; NaN];
        pos = pos + 1;
        err(1) = true;
        continue
    end
    if ~isempty(homologs)
        SA(homologs, 3) = 0;  % mask by 0 suffix frequency
        filtered = filtered + 1;
    end

    [B{blk, 3}, B{blk, 1}, B{blk, 2}] = most_freq_prefix(blockAA, SA, refAA, refNT);
    cmap_origin(pos:pos+m-1, 1) = B{blk, 1};
    cmap_origin(pos:pos+m-1, 2) = blk;

    same = cmap_origin(:, 1) == cmap_origin(pos, 1);
    if (mean(same) > max_pos) && (n > 1)
        SA(SA(:, 2) == cmap_origin(pos, 1), 3) = 0;  % mask
        filtered = filtered + 1;
        pos = find(same, 1);  % backtrack
        blk = cmap_origin(pos, 2) - 1;
        [B{blk+1:end, :}] = deal({});
        cmap_origin(pos:end, :) = 0;
    else
        pos = pos + m;
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
