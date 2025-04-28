#!/bin/bash

#Clear output
rm -f *output.txt

#Subsetting based on landcodes
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
#VCF file 
VCF="$HOME/Desktop/Stockholm_things/Basic_data/AADR/v62.0_1240k_public.vcf.gz"

#Output location 
Output_location="$HOME/Desktop/Stockholm_things/FUT2/Code/"

#Variants+Proxies used 
Variant="19:49206674-49206674"
Proxy1="19:49206417-49206417"
Proxy2="19:49206172-49206172"

#Loop
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

	#Write header
	echo -e "START\tSTOP\tN_IND\tN_CHR_target\tAF_target\tN_CHR_proxy1\tAF_proxy1\tN_CHR_proxy2\tAF_proxy2\tN_CHR_combined\tAF_combined" >> ${j}_output.txt
	tail -n 1 ${j}_output.txt
	rm -f samples1.txt		

	#Present day, $10 is time period
	awk  -F "\t" '{if ($10==0 && '"$location"') print $1}' "$M49Anno" >> samples1.txt

	x=($(wc samples1.txt))
	n_target=${x[0]}
	if [ ! -s samples1.txt ] 
			#The file is empty 
		then
			echo -e "0\t0\t$n_target\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan" >> ${j}_output.txt
		else 
			#The file is not empty
			#Create variables of values of the number of chromosomes ($4) and allele frequency ($8)
			n_CHR_target=($(bcftools view $VCF -S samples1.txt $Variant | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $4/2}'))
			AF_target=($(bcftools view $VCF -S samples1.txt $Variant | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $8}'))
	
				#Check first proxy
				bcftools view $VCF -S samples1.txt $Variant  2> /dev/null | tail -n 2 | awk '{for (i=1;i<=NF;i++) printf $i "\t"}' | awk '{for (i=1;i<=NF;i++) if ($i=="./.") print $(i-NF/2)}' > samples1_missing_gt_1.txt
				x=($(wc samples1_missing_gt_1.txt))
				n_proxy1=${x[0]}

				if [ ! -s samples1_missing_gt_1.txt ] 
					# The file is empty, and hence there are no individuals after subsetting
					then
						let n_CHR_proxy1=0 
						AF_proxy1="nan" 
					else
						n_CHR_proxy1=($(bcftools view $VCF -S samples1_missing_gt_1.txt $Proxy1  2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $4/2}'))
						AF_proxy1=($(bcftools view $VCF -S samples1_missing_gt_1.txt $Proxy1  2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $8}'))
				fi		
		
				#Check second proxy 
			bcftools view $VCF -S samples1_missing_gt_1.txt $Proxy1 2> /dev/null | tail -n 2 | awk '{for (i=1;i<=NF;i++) printf $i "\t"}' | awk '{for (i=1;i<=NF;i++) if ($i=="./.") print $(i-NF/2)}' > samples1_missing_gt_2.txt
			x=($(wc samples1_missing_gt_2.txt))
			n_proxy2=${x[0]}

			if [ ! -s samples1_missing_gt_2.txt ] 
				then
					let n_CHR_proxy2=0
					AF_proxy2="nan" 
				else
					n_CHR_proxy2=($(bcftools view $VCF -S samples1_missing_gt_2.txt $Proxy2  2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $4/2}'))
					AF_proxy2=($(bcftools view $VCF -S samples1_missing_gt_2.txt $Proxy2  2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $8}'))
			fi

			# Generate total number of chromosomes
			n_CHR_combined=$(( $n_CHR_target+$n_CHR_proxy1+$n_CHR_proxy2 ))

			if (( $n_CHR_combined >0 ))   
				then 
				#Generate the combined allele frequency of proxies and target variants with weighting based on number of chromosomes 
				AF_combined=$(echo "scale=5; ($n_CHR_target*$AF_target+$n_CHR_proxy1*$AF_proxy1+$n_CHR_proxy2*$AF_proxy2)/$n_CHR_combined" | bc)
				echo -e "0\t0\t$n_target\t$n_CHR_target\t$AF_target\t$n_CHR_proxy1\t$AF_proxy1\t$n_CHR_proxy2\t$AF_proxy2\t$n_CHR_combined\t$AF_combined" >> ${j}_output.txt
				else 
				#If total number of chromosomes was 0, would generate error message as division by 0
				echo -e "0\t0\t$n_target\t$n_CHR_target\t$AF_target\t$n_CHR_proxy1\t$AF_proxy1\t$n_CHR_proxy2\t$AF_proxy2\t$n_CHR_combined\tnan" >> ${j}_output.txt
			fi	
			tail -n 1 ${j}_output.txt
	fi


	#Ancient in windows	
	for i in {1..1000}
	do
		#Extract out samples in windows of 1000 years ($10 is the time period) in the areas we are interested in
		rm -f samples.txt	
		awk -F "\t" -v var="$i" '{if ($10>(10*var-500) && $10<=(10*var+500) && '"$location"') print $1}' "$M49Anno" > samples.txt 
		x=($(wc samples.txt))
		n_target=${x[0]}
		start_year=$(( 10*($i)-500 ))
		end_year=$(( 10*($i)+500 ))

		if [ ! -s samples.txt ]
			then
				# The file is empty with no individuals in this time period and area
				echo -e "$start_year\t$end_year\t$n_target\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan" >> ${j}_output.txt

			else 
				# The file is not empty
				# create variables of values of the number of chromosomes ($4) ($4/2 as pseudohaploid data) and allele frequency ($8)
				n_CHR_target=($(bcftools view $VCF -S samples.txt $Variant  | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $4/2}'))
				AF_target=($(bcftools view $VCF -S samples.txt $Variant | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $8}'))
	
				#Check first proxy
				bcftools view $VCF -S samples.txt $Variant  2> /dev/null | tail -n 2 | awk '{for (i=1;i<=NF;i++) printf $i "\t"}' | awk '{for (i=1;i<=NF;i++) if ($i=="./.") print $(i-NF/2)}' > samples_missing_gt_1.txt
				x=($(wc samples_missing_gt_1.txt))
				n_proxy1=${x[0]}

				if [ ! -s samples_missing_gt_1.txt ] 
					# The file is empty, and hence there are no individuals after subsetting
					then
						let n_CHR_proxy1=0 
						AF_proxy1="nan" 
					else
						n_CHR_proxy1=($(bcftools view $VCF -S samples_missing_gt_1.txt $Proxy1  2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $4/2}'))
						AF_proxy1=($(bcftools view $VCF -S samples_missing_gt_1.txt $Proxy1  2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $8}'))
				fi		
		
				#Check second proxy 
				bcftools view $VCF -S samples_missing_gt_1.txt $Proxy1 2> /dev/null | tail -n 2 | awk '{for (i=1;i<=NF;i++) printf $i "\t"}' | awk '{for (i=1;i<=NF;i++) if ($i=="./.") print $(i-NF/2)}' > samples_missing_gt_2.txt
				x=($(wc samples_missing_gt_2.txt))
				n_proxy2=${x[0]}

				if [ ! -s samples_missing_gt_2.txt ] 
					then
						let n_CHR_proxy2=0
						AF_proxy2="nan" 
					else
						n_CHR_proxy2=($(bcftools view $VCF -S samples_missing_gt_2.txt $Proxy2 2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $4/2}'))
						AF_proxy2=($(bcftools view $VCF -S samples_missing_gt_2.txt $Proxy2  2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $8}'))
				fi

				n_CHR_combined=$(( $n_CHR_target+$n_CHR_proxy1+$n_CHR_proxy2 )) 

				if (( $n_CHR_combined >0 ))   
					then 
					AF_combined=$(echo "scale=5; ($n_CHR_target*$AF_target+$n_CHR_proxy1*$AF_proxy1+$n_CHR_proxy2*$AF_proxy2)/$n_CHR_combined" | bc)
					echo -e "$start_year\t$end_year\t$n_target\t$n_CHR_target\t$AF_target\t$n_CHR_proxy1\t$AF_proxy1\t$n_CHR_proxy2\t$AF_proxy2\t$n_CHR_combined\t$AF_combined" >> ${j}_output.txt
					else 
					echo -e "$start_year\t$end_year\t$n_target\t$n_CHR_target\t$AF_target\t$n_CHR_proxy1\t$AF_proxy1\t$n_CHR_proxy2\t$AF_proxy2\t$n_CHR_combined\tnan" >> ${j}_output.txt
				fi
	
		fi
		tail -n 1 ${j}_output.txt

	done

	rm -f samples2.txt
	#Most ancient DNA >10000 years
	awk -F "\t" '{if ($10>10000 && '"$location"') print $1}' "$M49Anno" >> samples2.txt

	x=($(wc samples2.txt))
	n_target=${x[0]}
	if [ ! -s samples2.txt ] 
		then
			echo -e "10000\t50000\t$n_target\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan" >> ${j}_output.txt
		else 
			n_CHR_target=($(bcftools view $VCF -S samples2.txt $Variant  | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $4/2}'))
			AF_target=($(bcftools view $VCF -S samples2.txt $Variant | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $8}'))
	
				#Check first proxy
				bcftools view $VCF -S samples2.txt $Variant  2> /dev/null | tail -n 2 | awk '{for (i=1;i<=NF;i++) printf $i "\t"}' | awk '{for (i=1;i<=NF;i++) if ($i=="./.") print $(i-NF/2)}' > samples2_missing_gt_1.txt
				x=($(wc samples2_missing_gt_1.txt))
				n_proxy1=${x[0]}

				if [ ! -s samples2_missing_gt_1.txt ] 
					# The file is empty, and hence there are no individuals after subsetting
					then
						let n_CHR_proxy1=0 
						AF_proxy1="nan" 
					else
						n_CHR_proxy1=($(bcftools view $VCF -S samples2_missing_gt_1.txt $Proxy1 2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $4/2}'))
						AF_proxy1=($(bcftools view $VCF -S samples2_missing_gt_1.txt $Proxy1  2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $8}'))
				fi		
		
				#Check second proxy 
				bcftools view $VCF -S samples2_missing_gt_1.txt $Proxy1 2> /dev/null | tail -n 2 | awk '{for (i=1;i<=NF;i++) printf $i "\t"}' | awk '{for (i=1;i<=NF;i++) if ($i=="./.") print $(i-NF/2)}' > samples2_missing_gt_2.txt
				x=($(wc samples2_missing_gt_2.txt))
				n_proxy2=${x[0]}

				if [ ! -s samples2_missing_gt_2.txt ] 
					then
						let n_CHR_proxy2=0
						AF_proxy2="nan" 
					else
						n_CHR_proxy2=($(bcftools view $VCF -S samples2_missing_gt_2.txt $Proxy2  2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $4/2}'))
						AF_proxy2=($(bcftools view $VCF -S samples2_missing_gt_2.txt $Proxy2  2> /dev/null | vcftools --vcf - --freq --stdout | sed $'s/:/\t/g'  | tail -n +2 | awk '{print $8}'))
				fi

			n_CHR_combined=$(( $n_CHR_target+$n_CHR_proxy1+$n_CHR_proxy2 )) 

			if (( $n_CHR_combined >0 ))   
				then 
				AF_combined=$(echo "scale=5; ($n_CHR_target*$AF_target+$n_CHR_proxy1*$AF_proxy1+$n_CHR_proxy2*$AF_proxy2)/$n_CHR_combined" | bc)
				echo -e "10000\t50000\t$n_target\t$n_CHR_target\t$AF_target\t$n_CHR_proxy1\t$AF_proxy1\t$n_CHR_proxy2\t$AF_proxy2\t$n_CHR_combined\t$AF_combined" >> ${j}_output.txt
				else 
				echo -e "10000\t50000\t$n_target\t$n_CHR_target\t$AF_target\t$n_CHR_proxy1\t$AF_proxy1\t$n_CHR_proxy2\t$AF_proxy2\t$n_CHR_combined\tnan" >> ${j}_output.txt
			fi
	fi
	tail -n 1 ${j}_output.txt
	cp ${j}_output.txt $Output_location
done

