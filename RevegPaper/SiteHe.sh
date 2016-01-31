# Calculate Population Expected Heterozygosity from, per locus, allele frequencies (2 x ref allele x alt allele)
# Input: individual population .frq file (from vcftools)
# Output: '.frq_prefix'.He = modified .frq file to include additional column at end of (per locus) expected heterozygosity
# Usage: SiteHe.sh <input.frq>

awk 'NR==1 {print $0 "\t\t" "He"}
{if (NR >= 2) 
{split($5,ref,":"); split($6,alt,":"); print $0,2*ref[2]*alt[2]}
}' $1 > ${1%.frq}.He

