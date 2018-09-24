function codonseq = nt2codon(seq)
% codonseq = NT2CODON(seq)
% converts a sequence in NT alphabet and returns a sequence in codon
% alphabet, which is defined as follows:
% the index of the character at a position equals the index of the codon in
% the sorted list of all 64 triplets.
% this allows for all string algorithms, including the Chimera approach, to
% work seemlessly.
% ignores partial codons.
%
% Alon Diament, Tuller Lab, August 2018.

if iscellstr(seq)
    codonseq = cellfun(@nt2codon, seq, 'UniformOutput', false);
    return
end

if isempty(seq)
    codonseq = '';
    return
end

alphabet = fieldnames(codoncount(''));

n = length(seq);
n = n - mod(n, 3);
seq = cellstr(reshape(seq(1:n), 3, [])');

[~, ind] = ismember(seq, alphabet);

codonseq = char(ind)';
