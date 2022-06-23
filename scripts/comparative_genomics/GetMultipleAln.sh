Ref=[Referece genome]
CPU=[Number]
while read TE ; do   
    echo "Now looking for " $TE   
    sID=$(basename $TE | sed 's/#Unknown//')   
    grep -w "${TE}"  [bed file positions TEs].bed > $sID.gff
    cut -f 1,4-5 $sID.gff > $sID.bed
    python AddNToBed.py --N 2000 --b ${sID}.bed --o ${sID}_2000.bed 
    bedtools getfasta -fi $Ref -bed ${sID}_2000.bed -fo $sID.fa
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  $sID.fa |  awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' | sort -k1r,1n  | cut -f 2- | tr "\t" "\n" > $sID.sorted.fa
    head -n 200 $sID.sorted.fa > $sID.Long.fa
    clustalo -i $sID.Long.fa -o $sID.clo -t DNA --outfmt=clu --full --dealign --threads=$CPU
    echo "Done with " $sID
done < TE.fofn