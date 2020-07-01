#!/bin/bash
#
# #adapted from http://merenlab.org/2016/06/22/anvio-tutorial-v2
#
# # path="/media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/kbase_assembly_outputs/merged_MEGAHIT_assembly/"
# output="/media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/anvio_sgg_mag_processing/"
# reads="/media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/qc_trimmed_data/"
#
# #re-formatting input fastas
# # for spl in `cat samples.txt`; do
# 	# spl=$(echo $dir | sed -e "s/.FASTA//g")
# 	# spl_path=$(echo $path"$dir""/""$spl"".fa")
# 	# echo $spl
# 	# echo $output"$spl""-fixed.fa"
# 	# echo $output"$spl"".db"
#     #     echo
#     #     echo "Re-formatting "$spl
#     #     echo
#  	# 	anvi-script-reformat-fasta \
# 		# $spl_path \
# 		# -o $output"$spl""-fixed.fa" \
# 		# -l 0 --simplify-names
# echo
# echo "Creating contigs database"
# echo
# anvi-gen-contigs-database \
# 	-f final_contigs-fixed.fa_assembly.fa \
# 	-o $output"contigs.db" \
# 	-n "SGG metagenomes database"
# echo
# echo "Running HHMs"
# echo
# anvi-run-hmms \
# 	-c $output"contigs.db" \
#   --num-threads 8
# echo
# echo "Getting contig stats"
# echo
# anvi-script-gen_stats_for_single_copy_genes.py \
# 	$output"contigs.db"
# anvi-script-gen_stats_for_single_copy_genes.R \
# 	$output"contigs.db.hits" \
# 	$output"contigs.db.genes"
# echo
# echo "Setting up NCBI COG annotation thingy"
# echo
# anvi-setup-ncbi-cogs \
#   -T 8 --just-do-it
# echo
# echo "Adding NCBI COGS to contigs.db"
# anvi-run-ncbi-cogs \
#   -c $output"contigs.db" \
#   -T 8
# echo
# echo "Adding taxonomy with kaiju"
# echo
# anvi-get-sequences-for-gene-calls \
#   -c $output"contigs.db" \
#   -o $output"contigs_gene_calls.fa"
#
# kaiju -t /media/andre/B2F8C9A0F8C962E9/kaiju_proGenomes_db/nodes.dmp \
#   -f /media/andre/B2F8C9A0F8C962E9/kaiju_proGenomes_db/kaiju_db.fmi \
#   -i $output"contigs_gene_calls.fa" \
#   -o $output"contigs_gene_calls_nr_out" \
#   -z 8 \
#   -v
#
# addTaxonNames -t /media/andre/B2F8C9A0F8C962E9/kaiju_proGenomes_db/nodes.dmp \
#   -n /media/andre/B2F8C9A0F8C962E9/kaiju_proGenomes_db/names.dmp \
#   -i $output"contigs_gene_calls_nr_out" \
#   -o $output"contigs_gene_calls_nr_out.names" \
#   -r superkingdom,phylum,order,class,family,genus,species
#
# anvi-import-taxonomy-for-genes \
#   -i $output"contigs_gene_calls_nr_out.names" \
#   -c $output"contigs.db" \
#   -p kaiju \
# 	--just-do-it
# echo
# echo "Creating bowtie2 index"
# echo
# bowtie2-build \
# 	final_contigs-fixed.fa_assembly.fa \
# 	final_contigs-fixed.fa_assembly
#
# for spl in `cat samples.txt`; do
# 	echo
# 	echo "Adding BAM files to contigs database"
# 	echo
# 	echo "bowtie2 running for "$spl
# 	# gunzip -k $R1R2"$spl"_R1_001.fastq.gz $R1R2"$spl"_R1_001.fastq
# 	bowtie2 -x final_contigs-fixed.fa_assembly -q \
# 		-1 $reads"$spl"_trimqc_R1_001.fastq.gz \
# 		-2 $reads"$spl"_trimqc_R2_001.fastq.gz \
# 		--no-unal -p 8 -S $spl".sam"
# 	# echo $reads"$spl"_R1_001.fastq.gz
# 	echo
# 	echo "samtools view running for "$spl
# 	echo
# 	samtools view -bS -o $spl"-raw.bam" \
# 		$spl".sam"
# 	if [ ! -f $spl"-raw.s.bam" ]; then
# 		echo
# 		echo "samtools sort running for "$spl
# 		echo
#     samtools sort $spl"-raw.bam" -o $spl"-raw.s.bam" > $spl"-raw.s.bam"
# 		echo
# 		echo "samtools index running for "$spl
# 		echo
# 		samtools index $spl"-raw.s.bam"
# 	fi
# 	echo
# 	echo "anvi-profile for "$spl
# 	echo
# 	#not clustering contigs
#   anvi-profile \
#     -i $spl"-raw.s.bam" \
#     -c $output"contigs.db" \
# 		--output-dir $spl \
# 		--sample-name $spl \
# 		--min-contig-length 1000 \
# 		-W -T 8
# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest \
# 	anvi-migrate-db	$spl/PROFILE.db
# done
# # echo
# # echo "Merging all contigs databases"
# # echo
# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest \
# 	anvi-merge \
# 		S1/PROFILE.db S2/PROFILE.db S3/PROFILE.db \
# 		S4/PROFILE.db S5/PROFILE.db S6/PROFILE.db \
# 		S7/PROFILE.db S8/PROFILE.db S9/PROFILE.db \
# 		-o just_to_check_SAMPLES-MERGED \
# 		-c $output"contigs.db"
# echo
# echo "Generating summary files for SGG"
# echo
# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest \
# 	anvi-summarize \
# 		-p just_to_check_SAMPLES-MERGED/PROFILE.db \
# 		-c $output"contigs.db" \
# 		-C "CONCOCT" \
# 		-o just_to_check_SAMPLES-MERGED_CONCOCT
# echo
# echo "Refining MAGs with anvi-refine"
# echo
# anvi-refine -p SAMPLES-MERGED/PROFILE.db -c contigs.db -C CONCOCT -b Bin_XX --taxonomic-level t_class

#####Adding Maxbin2 binning to the profile
# #export split contigs .fa to bin with MaxBin2 and import as another binning collection
# anvi-export-splits-and-coverages -c contigs.db \
		# -p SAMPLES-MERGED/PROFILE.db \
		# --splits-mode \
		# -o ./split_contigs_table
#kbase churning
# sh anvio_bin_import.sh
#jfc had to grab an anvi'o image from somewhere to keep running this
# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` meren/anvio:latest anvi-import-collection \
# 	kbase_MaxBin2_outputs__bin_list.tsv \
# 	-p SAMPLES-MERGED/PROFILE.db \
# 	-c contigs.db \
# 	-C "Maxbin2_refined"
# #check out Maxbin2 binning
# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` meren/anvio:latest anvi-summarize \
# 	-p /media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/anvio_sgg_mag_processing/SAMPLES-MERGED/PROFILE.db \
# 	-c /media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/anvio_sgg_mag_processing/contigs.db \
# 	-C "Maxbin2_refined" \
# 	-o /media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/anvio_sgg_mag_processing/MERGED_SUMMARY_Maxbin2_refined
	#same for concoct binning
	# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` meren/anvio:latest anvi-summarize \
	# 	-p /media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/anvio_sgg_mag_processing/SAMPLES-MERGED/PROFILE.db \
	# 	-c /media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/anvio_sgg_mag_processing/contigs.db \
	# 	-C "CONCOCT" \
	# 	-o /media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/anvio_sgg_mag_processing/MERGED_SUMMARY_CONCOCT_refined
#refining with docker image
# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest anvi-refine -p SAMPLES-MERGED/PROFILE.db -c contigs.db -C Maxbin2_refined -b Bin_056 --taxonomic-level t_genus

# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest \
# 	anvi-interactive \
# 		-p SAMPLES-MERGED/PROFILE.db \
# 		-c contigs.db \
# 		-C Maxbin2_refined

sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest \
	anvi-search-functions \
		-c contigs.db \
		--search-terms methan \
		--verbose

# #filtering genomes for phylogenomics within anvi'o
# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest anvi-get-sequences-for-hmm-hits \
# 	-c /media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/anvio_sgg_mag_processing/contigs.db \
# 	-p /media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/anvio_sgg_mag_processing/SAMPLES-MERGED/PROFILE.db \
# 	--hmm-sources Campbell_et_al -o Campbell_et_al_seqs-for-phylogenomics.fa \
# 	--concatenate-genes --return-best-hit --get-aa-sequences \
# 	--gene-names Arg_tRNA_synt_N,B5,CTP_synth_N,CoaE,Competence,Cons_hypoth95,Cytidylate_kin,DNA_pol3_beta,DNA_pol3_beta_2,DNA_pol3_beta_3,EF_TS,Enolase_C,Enolase_N,FAD_syn,FDX-ACB,Flavokinase,GAD,GMP_synt_C,GTP1_OBG,GidB,GrpE,IF-2,IF2_N,IF3_C,IF3_N,IPPT,LepA_C,Methyltransf_5,MurB_C,NusA_N,Oligomerisation,PGK,PNPase,Pept_tRNA_hydro,Peptidase_A8,Phe_tRNA-synt_N,PseudoU_synth_1,RBFA,RNA_pol_A_CTD,RNA_pol_A_bac,RNA_pol_L,RNA_pol_Rpb1_1,RNA_pol_Rpb1_2,RNA_pol_Rpb1_3,RNA_pol_Rpb1_4,RNA_pol_Rpb1_5,RNA_pol_Rpb2_1,RNA_pol_Rpb2_2,RNA_pol_Rpb2_3,RNA_pol_Rpb2_45,RNA_pol_Rpb2_6,RNA_pol_Rpb2_7,RRF,RecA,RecR,Ribonuclease_P,Ribosom_S12_S23,Ribosomal_L1,Ribosomal_L10,Ribosomal_L11,Ribosomal_L11_N,Ribosomal_L12,Ribosomal_L13,Ribosomal_L14,Ribosomal_L16,Ribosomal_L17,Ribosomal_L18e,Ribosomal_L18p,Ribosomal_L19,Ribosomal_L2,Ribosomal_L20,Ribosomal_L21p,Ribosomal_L22,Ribosomal_L23,Ribosomal_L27,Ribosomal_L28,Ribosomal_L29,Ribosomal_L2_C,Ribosomal_L3,Ribosomal_L32p,Ribosomal_L35p,Ribosomal_L4,Ribosomal_L5,Ribosomal_L5_C,Ribosomal_L6,Ribosomal_L9_C,Ribosomal_L9_N,Ribosomal_S10,Ribosomal_S11,Ribosomal_S13,Ribosomal_S15,Ribosomal_S16,Ribosomal_S17,Ribosomal_S18,Ribosomal_S19,Ribosomal_S2,Ribosomal_S20p,Ribosomal_S3_C,Ribosomal_S4,Ribosomal_S5,Ribosomal_S5_C,Ribosomal_S6,Ribosomal_S7,Ribosomal_S8,Ribosomal_S9,RimM,RuvA_C,RuvA_N,RuvB_C,S-AdoMet_synt_C,S-AdoMet_synt_M,SRP_SPB,SecE,SecG,SecY,Seryl_tRNA_N,SmpB,THF_DHG_CYH,THF_DHG_CYH_C,TIM,TRCF,Toprim_N,Trigger_C,Trigger_N,TruB_N,UBA,UPF0054,UPF0079,UPF0081,UvrB,UvrC_HhH_N,Val_tRNA-synt_C,YchF-GTPase_C,dsrm,eIF-1a,tRNA-synt_1d,tRNA-synt_2d,tRNA_m1G_MT,zf-C4_ClpX
#
# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest anvi-gen-phylogenomic-tree \
# 	-f Campbell_et_al_seqs-for-phylogenomics.fa \
# 	-o Campbell_et_al_phylogenomics_tree.txt

	# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest anvi-get-sequences-for-hmm-hits \
	# 	-i internal_genomes.txt \
	# 	-p /media/andre/B2F8C9A0F8C962E9/SGG_metagenomes/data/anvio_sgg_mag_processing/SAMPLES-MERGED/PROFILE.db \
	# 	--hmm-sources Rinke_et_al -o Rinke_et_al_seqs-for-phylogenomics.fa \
	# 	--concatenate-genes --return-best-hit --get-aa-sequences \
	# 	--gene-names ATP-synt_C,ATP-synt_D,ATP-synt_F,Adenylsucc_synt,AdoHcyase,AdoHcyase_NAD,AdoMet_Synthase,Archease,B5,CTP-dep_RFKase,CTP_synth_N,DALR_1,DFP,DHO_dh,DKCLD,DNA_binding_1,DNA_primase_S,DNA_primase_lrg,DUF137,DUF357,DUF359,DUF46,DUF655,DUF814,DUF99,Diphthamide_syn,EF1_GNE,EFG_C,EFG_IV,EIF_2_alpha,Enolase_C,Fibrillarin,GAD,GCD14,GMP_synt_C,HMG-CoA_red,Ham1p_like,IF-2,LigT_PEase,MAF_flag10,MoaC,Mob_synth_C,NAC,NDK,NMD3,Nop,Nop10p,PGK,PTH2,PcrB,Plug_translocon,Prefoldin,Prefoldin_2,PyrI,PyrI_C,RNA_pol_A_bac,RNA_pol_N,RNA_pol_Rpb1_1,RNA_pol_Rpb1_2,RNA_pol_Rpb1_3,RNA_pol_Rpb1_4,RNA_pol_Rpb2_1,RNA_pol_Rpb2_2,RNA_pol_Rpb2_3,RNA_pol_Rpb2_4,RNA_pol_Rpb2_5,RNA_pol_Rpb2_6,RNA_pol_Rpb2_7,RNA_pol_Rpb4,RNA_pol_Rpb5_C,RNA_pol_Rpb6,RNase_HII,RS4NT,Rib_5-P_isom_A,Ribosom_S12_S23,Ribosomal_L1,Ribosomal_L10,Ribosomal_L11,Ribosomal_L11_N,Ribosomal_L13,Ribosomal_L14,Ribosomal_L15e,Ribosomal_L16,Ribosomal_L18p,Ribosomal_L19e,Ribosomal_L2,Ribosomal_L21e,Ribosomal_L22,Ribosomal_L23,Ribosomal_L24e,Ribosomal_L29,Ribosomal_L2_C,Ribosomal_L3,Ribosomal_L30,Ribosomal_L31e,Ribosomal_L32e,Ribosomal_L37ae,Ribosomal_L37e,Ribosomal_L39,Ribosomal_L4,Ribosomal_L44,Ribosomal_L5,Ribosomal_L5_C,Ribosomal_L6,Ribosomal_S11,Ribosomal_S13,Ribosomal_S13_N,Ribosomal_S14,Ribosomal_S15,Ribosomal_S17,Ribosomal_S17e,Ribosomal_S19,Ribosomal_S19e,Ribosomal_S2,Ribosomal_S24e,Ribosomal_S27,Ribosomal_S27e,Ribosomal_S28e,Ribosomal_S3Ae,Ribosomal_S3_C,Ribosomal_S4e,Ribosomal_S5,Ribosomal_S5_C,Ribosomal_S6e,Ribosomal_S7,Ribosomal_S8,Ribosomal_S8e,Ribosomal_S9,RtcB,SBDS,SBDS_C,SHS2_Rpb7-N,SNO,SRP19,SRP_SPB,SUI1,Sec61_beta,SecY,Spt4,Spt5-NGN,TIM,TP6A_N,TRM,Topo-VIb_trans,TruB_N,UPF0004,UPF0086,V_ATPase_I,Wyosine_form,XPG_I,XPG_N,YjeF_N,dsDNA_bind,eIF-5a,eIF-6,eIF2_C,tRNA-synt_1c,tRNA-synt_1c_C,tRNA-synt_1d,tRNA_NucTransf2,tRNA_deacylase,vATP-synt_E
	#
	# sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest anvi-gen-phylogenomic-tree \
	# 	-f Rinke_et_al_seqs-for-phylogenomics.fa \
	# 	-o Rinke_et_al_phylogenomics_tree.txt
