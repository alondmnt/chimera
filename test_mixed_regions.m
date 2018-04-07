chim = fastaread('../../output-chim.fasta');
cod = fastaread('../../output-cod.fasta');
def = fastaread('../../output-def.fasta');

mix = fastaread('../../output-mix.fasta');
mix_s = fastaread('../../output-mix-single.fasta');

assert(strcmp(mix(2).Sequence, mix_s.Sequence));

chim_pos = cellfun(@(x, y) {x(1:3:end) == y(1:3:end)}, {mix.Sequence}, {chim.Sequence});
cod_pos = cellfun(@(x, y) {x(1:3:end) == y(1:3:end)}, {mix.Sequence}, {cod.Sequence});
def_pos = cellfun(@(x, y) {x(1:3:end) == y(1:3:end)}, {mix.Sequence}, {def.Sequence});
