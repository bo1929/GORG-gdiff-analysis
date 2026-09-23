Sketch selected genomes:
```
/bin/gdiff-x86 sketch -l 0 -k 23 -w 23 --frac 0.5 --input-list selected_genomes-ref.txt -o selected_genomes.gdsk
```

Roll windows over selected:
```
cat selected_genomes-queries.txt  | cut -f2 -d'/' | sed 's/_contigs.fasta//' | xargs -P 16 -I{} bash -c "./bin/gdiff roll contigs-gt80-complete/{}_contigs.fasta selected_genomes.gdsk -l 500 -s 500  -o results/gdiff-roll/{}.tsv"
```
