makeblastdb -in all_plasmids.fasta -out blastout/all_plasmids -dbtype nucl
mkdir -p blastout
for f in assm/*.fasta 
do
	name=$(basename $f .fasta)
	blastn -query $f -db blastout/all_plasmids -outfmt '6 std qcovs slen' > blastout/$name.blast
done

Rscript makeTable.r

