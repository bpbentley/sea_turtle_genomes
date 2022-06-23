
OutREblast=[output reciprocal Blast]

cat $OutREblast | sed '{s/,/\t/g}' | sed '{s/-/\t/g}' > $OutREblast_mod.tsv

awk '{if ($8>'$STARTregion' && $9<'$ENDregion') print}' $OutREblast_mod.tsv | cut -f 10 | sort |uniq -c | awk '{if ($1 != 1) print}' | awk -F ":" '{print $2}' > [file of duplications]

file=$OutREblast_mod.tsv
while read -r line; do
        echo $line
        grep $line $file
done < [file with the duplicated IDs]