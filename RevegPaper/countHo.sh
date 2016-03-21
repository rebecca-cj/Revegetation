## Count number of homozygous and heterozygous loci per individual (Observed Heterozygosity) ##
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
# Input: .012 vcftools output file (rows = individuals, columns = loci)
# Output: '.012_prefix'.het = text file with one row per individual; columns being sample = sample number (from .012 file), ref = count homozygous (reference allele), het = count heterozygotes, alt = count homozygous (alternate allele), missing = count missing genotypes, total = total number of genotypes (incl. missing), propH = proportion of heterozgotes (based on loci with genotypes)
# Usage: countHo.sh <input.012>
#####

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
