for file in `find $(pwd) -name '*.faa' -not -path '*/\.*'`; 
do binID="`basename $file | cut -f1 -d '.'`";
geneCount="`grep -o '>' $file | wc -l`";
echo "${binID}	${geneCount}" >> geneCounts.txt;
done

for file in `find $(pwd) -name '*.fasta' -not -path '*/\.*'`; 
do binID="`basename $file | cut -f1 -d '.'`";
nucleotideCount="`grep -v '>' $file | wc -c`";
echo "${binID}	${nucleotideCount}" >> nucleotideCounts.txt;
done

