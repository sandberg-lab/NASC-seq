experimentInfo = dict(
	readlength='150',
	paired=True,
	verbose=False,
	gtf='/mnt/kauffman/hendgert/genomes/hg38_ercc/Homo_sapiens.GRCh38.92.ERCC.chr.gtf',
	readcountdir='star_hg38_ercc_readcounts',
	output='star_hg38_ercc_outfiles',
	aligndir='star_hg38_ercc',
	gnv='/mnt/kauffman/hendgert/genomes/hg38_ercc',
	gtf='/mnt/kauffman/hendgert/genomes/hg38_ercc/Homo_sapiens.GRCh38.92.ERCC.chr.gtf',
	sjdb='/mnt/kauffman/hendgert/genomes/hg38_ercc/sjdbList.fromGTF.out.tab'
	)

cutoffs = dict(
	totalReads=1000000,
	assignedReads=500000,
	percentageAssigned=50,
	numberofHighExprGenes=None,
	posratioCutoff=0.4,
	cellnumberCutoff=5
	)

distributions = dict(
	rootDir='/mnt/kauffman/hendgert/Programs/NASCseq',
	trimgaloreDist='/mnt/kauffman/hendgert/Programs/trimgalore/TrimGalore-0.4.5/trim_galore',
	starDist='/mnt/kauffman/chrisz/programs/STAR/bin/Linux_x86_64/STAR',
	picardDist='/mnt/kauffman/hendgert/Programs/picard/picard.jar',
	javaDist='/mnt/kauffman/chrisz/programs/jdk1.8.0_161/bin/java'
	)