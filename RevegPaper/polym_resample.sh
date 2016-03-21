## Calculate percentage polymorphic loci for random subset of samples (n), from x number of replicates (with replacement between replicates) ##
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
# Input: "xxx.family" = text file of sample names for a single population in plink format i.e. 1st column = 'Family', 2nd column = 'sample name'
# Input: genotypes in PED format (plink);  script input = prefix for .ped and .map files
# Output: "output_prefix"_poly."number"n.reps = 
# Usage: polym_resample.sh <xxx.family> <plink_file_prefix> <num_of_reps_(x)> <samples_per_rep_(n)> <output_file_prefix>
#####

pop=${1%.family}
for num in `seq 1 $3`
do
if [ $num -le $3 ]; then
	shuf -n $4 $1 > $pop.$4.tmp.samples;
	printf "POP\t$pop\tREP\t$num\t" >> $5_rep.$4n.log
	awk '{printf $2 "\t"}END{printf "\n"}' $pop.$4.tmp.samples >> $5_rep.$4n.log
	plink --file $2 --keep $pop.$4.tmp.samples --out $pop_$4 --freq 1>>$5_rep.$4n.log;
	awk -v popname=$pop -v rep=$num 'NR>1 {if($5>0) poly++}END{print popname "\t" rep "\t" poly}' $pop_$4.frq >> $5_poly.$4n.reps;
	rm $pop.$4.tmp.samples $pop_$4.frq
fi
done
