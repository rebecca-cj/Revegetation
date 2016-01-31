## Count number of homozygous and heterozygous loci per individual (Observed Heterozygosity) ##
# Input: .012 vcftools output file (rows = individuals, columns = loci)
# Output: '.012_prefix'.het = text file with one row per individual; columns being sample = sample number (from .012 file), ref = count homozygous (reference allele), het = count heterozygotes, alt = count homozygous (alternate allele), missing = count missing genotypes, total = total number of genotypes (incl. missing), propH = proportion of heterozgotes (based on loci with genotypes)
# Usage: countHo.sh <input.012>

awk '
BEGIN{print "sample \t ref \t het \t alt \t missing \t total \t propH"}
{
ref=0; het=0; alt=0; missing=0; total=0; 
for(i=2; i <=NF; i++) 
	{if($i == 0) ref++; 
		else if($i==1) het++;
			else if($i==2) alt++;
				else if($i==-1) missing++};
total = (ref + het + alt + missing);
propH = (het/(total-missing));
printf("%d\t%d\t%d\t%d\t%d\t%d\t%f\n",$1,ref,het,alt,missing,total,propH)
}' $1 > ${1%.012}.het
