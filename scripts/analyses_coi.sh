

mothur

#create logfile with all manip that are going to be done
set.logfile(name=bioman_coi.log)

#analyse reads (length, ambiguities) and count number of reads per group
summary.seqs(fasta=bioman_coi.fasta, processors=12)
count.groups(group=bioman_coi.groups)

#dereplicate
unique.seqs(fasta=bioman_coi.fasta)
summary.seqs(fasta=bioman_coi.unique.fasta, name=bioman_coi.names, processors=12)

#align sequences against the 12S rRNA database
align.seqs(fasta=bioman_coi.unique.fasta, reference=/share/databases/BOLD/COI-BOLD-may2018_good-noDot_subset.align, processors=12, flip=T)
summary.seqs(fasta=bioman_coi.unique.align, name=bioman_coi.names, processors=12)

#remove sequences not covering the "coi" region of 12S
screen.seqs(fasta=bioman_coi.unique.align, name=bioman_coi.names, group=bioman_coi.groups, start=726, end=1055, minlength=313, maxlength=313, maxambig=0, processors=12)
summary.seqs(fasta=bioman_coi.unique.good.align, name=bioman_coi.good.names, processors=12)
count.groups(group=bioman_coi.good.groups)

# remove columns that contain gap characters
filter.seqs(fasta=bioman_coi.unique.good.align, vertical=T, processors=12)
unique.seqs(fasta=bioman_coi.unique.good.filter.fasta, name=bioman_coi.good.names)
summary.seqs(fasta=bioman_coi.unique.good.filter.unique.fasta,name=bioman_coi.unique.good.filter.names,processors=12)

#remove chimeras
chimera.uchime(fasta=bioman_coi.unique.good.filter.unique.fasta, name=bioman_coi.unique.good.filter.names, group=bioman_coi.good.groups, processors=12)
remove.seqs(accnos=bioman_coi.unique.good.filter.unique.denovo.uchime.accnos, fasta=bioman_coi.unique.good.filter.unique.fasta, name=bioman_coi.unique.good.filter.names, group=bioman_coi.good.groups, dups=T)
summary.seqs(fasta=bioman_coi.unique.good.filter.unique.pick.fasta, name=bioman_coi.unique.good.filter.pick.names, processors=12)
count.groups(group=bioman_coi.good.pick.groups)

#for clarity, rename files
system(cp bioman_coi.unique.good.filter.unique.pick.fasta bioman_coi_all.fasta)
system(cp bioman_coi.unique.good.filter.pick.names bioman_coi_all.names)
system(cp bioman_coi.good.pick.groups bioman_coi_all.groups)

#create count_table with the new files
make.table(name=bioman_coi_all.names, group=bioman_coi_all.groups, processors=12)

quit

# Clustering and taxonomic assignment

# VSEARCH + LULU
mothur "#cluster(fasta=bioman_coi_all.fasta, count=bioman_coi_all.count_table, method=agc)"
mv bioman_coi_all.agc.unique_list.list bioman_coi_all_vsearch.list
mothur "#get.oturep(list=bioman_coi_all_vsearch.list, fasta=bioman_coi_all.fasta, method=abundance, name=bioman_coi_all.names)"
awk '/^>/{printf ">Otu%05d\n",++i; next}{print}' bioman_coi_all_vsearch.0.03.rep.fasta | sed s/"\."//g| sed s/"-"//g > bioman_coi_all_vsearch.0.03.rep_otuName.fasta 
mothur "#make.shared(list=bioman_coi_all_vsearch.list, count=bioman_coi_all.count_table, label=0.03)"
makeblastdb -in bioman_coi_all_vsearch.0.03.rep_otuName.fasta -dbtype nucl
blastn -db bioman_coi_all_vsearch.0.03.rep_otuName.fasta -outfmt '6 qseqid sseqid pident' -out bioman_coi_all_vsearch.0.03.rep_otuName.dist -qcov_hsp_perc 100 -perc_identity 90 -query bioman_coi_all_vsearch.0.03.rep_otuName.fasta -num_threads 8
cut -f 2,4- bioman_coi_all_vsearch.shared > bioman_coi_all_vsearch_0.03.tsv
Rscript LULU2.R bioman_coi_all_vsearch.0.03.rep_otuName.dist bioman_coi_all_vsearch_0.03.tsv 97
head -1 bioman_coi_all_vsearch_0.03_curated.shared | cut -f 4- | sed s/"\t"/"\n"/g >  bioman_coi_all_vsearch_0.03_curated.outs
nO=(`wc bioman_coi_all_vsearch_0.03_curated.outs`)
n=`expr "$nO" + 3`
for f in `seq 4 $n`; do awk -v f="$f" '{s+=$f} END {print s}' bioman_coi_all_vsearch_0.03_curated.shared; done > bioman_coi_all_vsearch_0.03_curated.nRead_otu
grep -f bioman_coi_all_vsearch_0.03_curated.outs -A1 bioman_coi_all_vsearch.0.03.rep_otuName.fasta >  bioman_coi_all_vsearch.0.03.rep_otuName_curated.fasta
mothur "#classify.seqs(fasta=bioman_coi_all_vsearch.0.03.rep_otuName_curated.fasta, template=/share/databases/BOLD/COI-BOLD-may2018_good.fasta, taxonomy=/share/databases/BOLD/COI-BOLD-may2018_good.tax, method=wang, cutoff=30, processors=12)"
echo -e "OTU\tSize\tTaxonomy" > bioman_coi_all_vsearch.0.03.rep_otuName_curated.COI_BOLD_may2018_good.wang_corr.taxonomy
paste <(cut -f1 bioman_coi_all_vsearch.0.03.rep_otuName_curated.COI_BOLD_may2018_good.wang.taxonomy) <(cut -f1 bioman_coi_all_vsearch_0.03_curated.nRead_otu) <(cut -f 2- bioman_coi_all_vsearch.0.03.rep_otuName_curated.COI_BOLD_may2018_good.wang.taxonomy) | sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g >> bioman_coi_all_vsearch.0.03.rep_otuName_curated.COI_BOLD_may2018_good.wang_corr.taxonomy


# SWARM + LULU
Rscript PREP_SWARM.R bioman_coi_all.names bioman_coi_all.fasta
swarm bioman_coi_all_swarm.fasta -d 1 -t 10 -f -l bioman_coi_all_swarm_d1.log -r -o bioman_coi_all_swarm_d1.out -s bioman_coi_all_swarm_d1.stats -w bioman_coi_all_swarm_d1_rep.fasta 
awk '/^>/{printf ">Otu%04d\n",++i; next}{print}' bioman_coi_all_swarm_d1_rep.fasta > bioman_coi_all_swarm_d1_rep_otuName.fasta
sed -E  s/\(_[0-9]*,\)/","/g bioman_coi_all_swarm_d1.out   | sed -E s/\(_[0-9]*\\t\)/"\\t"/g | sed -E  s/\(_[0-9]*$\)//g | sed s/swarm/swarm_1/g >  bioman_coi_all_swarm_d1.list
mothur "#make.shared(list=bioman_coi_all_swarm_d1.list, count=bioman_coi_all.count_table, label=swarm_1)"
makeblastdb -in bioman_coi_all_swarm_d1_rep_otuName.fasta -dbtype nucl
blastn -db bioman_coi_all_swarm_d1_rep_otuName.fasta -outfmt '6 qseqid sseqid pident' -out bioman_coi_all_swarm_d1_rep_otuName.dist -qcov_hsp_perc 100 -perc_identity 90 -query bioman_coi_all_swarm_d1_rep_otuName.fasta -num_threads 8
cut -f 2,4- bioman_coi_all_swarm_d1.shared > bioman_coi_all_swarm_d1.tsv
Rscript LULU2.R bioman_coi_all_swarm_d1_rep_otuName.dist bioman_coi_all_swarm_d1.tsv 97
head -1 bioman_coi_all_swarm_d1_curated.shared | cut -f 4- | sed s/"\t"/"\n"/g >  bioman_coi_all_swarm_d1_curated.outs
nO=(`wc  bioman_coi_all_swarm_d1_curated.outs`)
n=`expr "$nO" + 3`
for f in `seq 4 $n`; do awk -v f="$f" '{s+=$f} END {print s}'  bioman_coi_all_swarm_d1_curated.shared; done >  bioman_coi_all_swarm_d1_curated.nRead_otu
grep -f bioman_coi_all_swarm_d1_curated.outs -A1 bioman_coi_all_swarm_d1_rep_otuName.fasta >  bioman_coi_swarm_d1_rep_otuName_curated.fasta
mothur "#classify.seqs(fasta=bioman_coi_swarm_d1_rep_otuName_curated.fasta, template=/share/databases/BOLD/COI-BOLD-may2018_good.fasta, taxonomy=/share/databases/BOLD/COI-BOLD-may2018_good.tax, method=wang, cutoff=30, processors=12)"
echo -e "OTU\tSize\tTaxonomy" > bioman_coi_swarm_d1_rep_otuName_curated.COI_BOLD_may2018_good.wang_corr.taxonomy
paste <(cut -f1 bioman_coi_swarm_d1_rep_otuName_curated.COI_BOLD_may2018_good.wang.taxonomy) <(cut -f1 bioman_coi_all_swarm_d1_curated.nRead_otu) <(cut -f 2- bioman_coi_swarm_d1_rep_otuName_curated.COI_BOLD_may2018_good.wang.taxonomy) | sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g >> bioman_coi_swarm_d1_rep_otuName_curated.COI_BOLD_may2018_good.wang_corr.taxonomy

	
# PHYLOTYPES (assign taxonomy without clustering)
mothur "#classify.seqs(fasta=bioman_coi_all.fasta, template=/share/databases/BOLD/COI-BOLD-may2018_good.fasta, taxonomy=/share/databases/BOLD/COI-BOLD-may2018_good.tax, group=bioman_coi_all.groups, method=wang, cutoff=25, processors=12)"
mothur "#phylotype(taxonomy=bioman_coi_all.COI_BOLD_may2018_good.wang.taxonomy)"
mothur "#make.shared(list=bioman_coi_all.COI_BOLD_may2018_good.wang.tx.list, count=bioman_coi_all.count_table, label=1)"
mothur "#classify.otu(list=bioman_coi_all.COI_BOLD_may2018_good.wang.tx.list, count=bioman_coi_all.count_table, taxonomy=bioman_coi_all.COI_BOLD_may2018_good.wang.taxonomy, label=1)"
sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g bioman_coi_all.COI_BOLD_may2018_good.wang.tx.1.cons.taxonomy> bioman_coi_all.COI_BOLD_may2018_good.wang.tx.1.cons_corr.taxonomy