## Create a list of loci that occur in at least one population (from the list given) with a minor allele frequency (maf) greater than given value ##
#
#  Copyright 2016 Rebecca Jordan
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#####
# Input = individual population .frq files (vcftools output) and text file containing population names eg. .frq file prefixes (minus .frq extension)
# Usage = filterSNPs.sh <pop.name.file>	<number.of.pops> <maf> <output.prefix>
#####

# create a list of biallelic loci found in in >80% of individuals in all populations (loci with greater than 2 alleles excluded)
	> $4_all.loci
	for pop in `cat $1`; do awk 'NR > 1 {if($3 <=2) print $1 "_" $2}' $pop".frq" >> $4_all.loci; done
	sort -n $4_all.loci | uniq -c - | awk -v number=$2 '{if($1 == number) print $0}' > $4_all.count

# filter each population for those loci found in all populations
	for pop in `cat $1`; do awk '{print $1 "_" $2 "\t" $0}' $pop".frq" | awk 'FNR==NR{a[$2]++;next}a[$1]' $4_all.count - > $4_${pop%.frq}.80loci; done

# reformat and score loci within each population
	for pop in `cat $1`; do 
		cut -f 2-7 $4_$pop".80loci" | 
		awk '{split($5,a,":"); split($6,b,":"); print $1 "_" $2 "\t" $3 "\t" $4 "\t" a[1] b[1] "\t" a[2] "\t" b[2]}' - | 
		awk -v maf=$3 '{if($6 == 1) print $1, $4, "alt"; 
					else if ($5 == 1) print $1, $4, "ref"; 
						else if ($5 >= maf && $6 >= maf) print $1, $4, "poly"; 
							else if($5 <= maf || $6 <= maf) print $1, $4, "na_MAF"
				}' - |
		awk 'BEGIN{print "Locus","Alleles","Pop"}{print}' - > $4_${pop%.80loci}.score
	done 
	for pop in `cat $1`; do sed -i "s/Pop/${pop%.score}/" $4_$pop".score"; done

# combine population score for loci
	head -1 $1 > tmp
	awk '{print $1}' $4_`cat tmp`".score" > $4_all.score
	for pop in `cat $1`; do
		awk '{print $1 "\t" $3}' $4_$pop".score" | join $4_all.score - > tmp.txt; mv tmp.txt $4_all.score; done
	rm tmp

# count scores per locus
	awk 'BEGIN{print "locus \t ref \t alt \t poly \t na_missing \t na_MAF \t total"}
	{ref=0; het=0; alt=0; poly=0; na_missing=0; na_MAF=0; total=0;
		for(i=2; i <=NF; i++) {
			if($i == "ref") ref++; 
				else if($i=="alt") alt++; 
				else if($i=="poly") poly++; 
				else if($i=="na_missing") na_missing++; 
				else if($i=="na_MAF") na_MAF++
			}; 
		total = (ref + alt + poly + na_missing + na_MAF); 
		printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\n",$1,ref,alt,poly,na_missing,na_MAF,total)}' $4_all.score > $4_all.MAFcount

# check for loci monomorphic within populations but potentially polymorphic between
	awk '{if($4==0 && $2>0 && $3 >0) print $0}' $4_all.MAFcount > $4_mono.loci

# make list of those loci polymorphic in at least one population
	awk 'BEGIN{NR = 1}{if($4>=1) print}' $4_all.MAFcount | awk 'NR > 1 {split($1,a,"_"); print a[1] "_" a[2] "\t" a[3]}' - > $4_keep_poly.loci
	
# make list of those loci polymorphic in at least one population AND all monomorphic loci
	awk 'BEGIN{NR = 1}{if($4>=1) print; else if($6==0 && $2 >=1) print; else if($6==0 && $3 >=1) print}' $4_all.MAFcount | awk 'NR > 1 {split($1,a,"_"); print a[1] "_" a[2] "\t" a[3]}' - > $4_keep_poly_mono.loci
	
# output counts
	wc -l $4_all.count $4_all.score $4_all.MAFcount $4_mono.loci $4_keep_poly.loci $4_keep_poly_mono.loci $4_*.80loci $4_*.score
