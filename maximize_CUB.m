function [seq_nt, CUB] = maximize_CUB(seq_aa, ref_nt)
% maximize_CUB(seq, ref)
%   replace all codons with optimal ones according to reference
%   sequences given in [ren_nt] (cell array). alternatively,
%   [ref_nt] can also be the struct output of codonbias().

% seq_aa = nt2aa(seq_nt, 'AlternativeStartCodons', false);
seq_nt = aa2nt(seq_aa);  % init

if isstruct(ref_nt)
    CUB = ref_nt;
else
%     test = ref_nt;
    for i = 1:length(ref_nt)
        last_legal_codon = length(ref_nt{i});
        last_legal_codon = last_legal_codon - mod(last_legal_codon, 3);
        ref_nt{i} = ref_nt{i}(1:last_legal_codon);  % length is now divisible by 3
    end
    CUB = codonbias(strcat(ref_nt{:}));
end

aa_list = fieldnames(CUB);
for a = 1:length(aa_list)
    aa = aminolookup(aa_list{a});
    pos = 3*(find(seq_aa == aa) - 1) + 1;  % aa to codon
    [~, bestcode] = max(CUB.(aa_list{a}).Freq);
    for p = pos
        seq_nt(p:p+2) = CUB.(aa_list{a}).Codon{bestcode};
    end
end

% test
if ~all(nt2aa(seq_nt, 'AlternativeStartCodons', false) == seq_aa)
    error('codon optimization error: too many cooks!');
end
