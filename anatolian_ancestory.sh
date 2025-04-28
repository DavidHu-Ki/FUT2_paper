#!/bin/bash
#subset based on M49 Anno file from AADR 
Britishisles="\$50==826"
Scandinavia="\$50==578 || \$50==752 || \$50==246 || \$50==208"
WestEu="\$50==250 || \$50==56 || \$50==528"
CentralEu="\$50==40 || \$50==203 || \$50==276 || \$50==348 || \$50==756 || \$50==616"
Iberian="\$50==724 || \$50==620"
Italy="\$50==380"
Balkan="\$50==100 || \$50==191 || \$50==300 || \$50==642 || \$50==688"
Europe="\$50==40 || \$50==203 || \$50==276 || \$50==348 || \$50==756 || \$50==616 || \$50==578 || \$50==752 || \$50==246 || \$50==208 || \$50==826 || \$50==250 || \$50==56 || \$50==528 || \$50==724 || \$50==620 || \$50==380 || \$50==100 || \$50==191 || \$50==300 || \$50==642 || \$50==688"
Africa="\$44==002"
Americas="\$44==019"
Asia="\$44==142"
EuropeandRussia="\$44==150"
Oceania="\$44==009"
EuropewestofUral="\$44==150 && \$18<=60"

#M49.Anno file 
M49Anno="$HOME/Desktop/Stockholm_things/Basic_data/AADR/v62.0_1240k_public_M49_het.anno"

Location_1="$Britishisles" 
Location_2="$Scandinavia"
Location_3="$WestEu"
Location_4="$CentralEu"
Location_5="$Iberian"
Location_6="$Italy"
Location_7="$Balkan"

Total_locations=7
for j in $(seq 1 $Total_locations); do
	location=$(eval echo "\$Location_$j") 
	echo "$location"
awk -F "\t" '{if ('"$location"') print $1}'"$M49Anno"   > ancestry_$j.txt 

grep -w -f "ancestry_$j.txt" Patterson_2022_ancestry.csv > ancestrycombined_$j.txt
done 
