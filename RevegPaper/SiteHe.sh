## Calculate Population Expected Heterozygosity from, per locus, allele frequencies (2 x ref allele x alt allele) ##
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
# Input: individual population .frq file (from vcftools)
# Output: '.frq_prefix'.He = modified .frq file to include additional column at end of (per locus) expected heterozygosity
# Usage: SiteHe.sh <input.frq>
#####

awk 'NR==1 {print $0 "\t\t" "He"}
{if (NR >= 2) 
{split($5,ref,":"); split($6,alt,":"); print $0,2*ref[2]*alt[2]}
}' $1 > ${1%.frq}.He

