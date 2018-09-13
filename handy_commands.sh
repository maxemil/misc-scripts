#Get all characters (eg check if sequences contain 'U')
grep -v '>' <fasta file> | od -cvAnone -w1 | sort -bu
