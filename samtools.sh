for f in S.sam
do
	echo $f
	samtools view -b $f > $f.bam
	echo 'view done, start sort'
	samtools sort $f.bam -o $f.sorted.bam
	echo 'sort done, start index'
	samtools index $f.sorted.bam
	echo 'index done, start idxstats'
	samtools idxstats $f.sorted.bam > $f.contigs.txt
done
