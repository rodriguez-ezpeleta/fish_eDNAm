
DBlocation=/share/projects/BIOMAN_eDNA/TELEO/local_DB/New_Database
DBprefix=rRNA_12S_teleo_worms_clean
globalDBlocation=/share/projects/BIOMAN_eDNA/TELEO/global_DB
globalDBprefix=global_teleo
mothur

#create logfile with all manip that are going to be done
set.logfile(name=bioman_teleo.log)

#analyse reads (length, ambiguities) and count number of reads per group
summary.seqs(fasta=bioman_teleo.fasta, processors=12)
count.groups(group=bioman_teleo.groups)

#dereplicate
unique.seqs(fasta=bioman_teleo.fasta)
summary.seqs(fasta=bioman_teleo.unique.fasta, name=bioman_teleo.names, processors=12)

#align sequences against the 12S rRNA database
align.seqs(fasta=bioman_teleo.unique.fasta, reference=../database/teleo.align, processors=12, flip=T)
summary.seqs(fasta=bioman_teleo.unique.align, name=bioman_teleo.names, processors=12)

#remove sequences not covering the "teleo" region of 12S
screen.seqs(fasta=bioman_teleo.unique.align, name=bioman_teleo.names, group=bioman_teleo.groups, start=6, end=77, minlength=60, maxlength=100, maxambig=0, processors=12)
summary.seqs(fasta=bioman_teleo.unique.good.align, name=bioman_teleo.good.names, processors=12)
count.groups(group=bioman_teleo.good.groups)

# remove columns that contain gap characters
filter.seqs(fasta=bioman_teleo.unique.good.align, vertical=T, processors=12)
unique.seqs(fasta=bioman_teleo.unique.good.filter.fasta, name=bioman_teleo.good.names)
summary.seqs(fasta=bioman_teleo.unique.good.filter.unique.fasta,name=bioman_teleo.unique.good.filter.names,processors=12)

#remove chimeras
chimera.uchime(fasta=bioman_teleo.unique.good.filter.unique.fasta, name=bioman_teleo.unique.good.filter.names, group=bioman_teleo.good.groups, processors=12)
remove.seqs(accnos=bioman_teleo.unique.good.filter.unique.denovo.uchime.accnos, fasta=bioman_teleo.unique.good.filter.unique.fasta, name=bioman_teleo.unique.good.filter.names, group=bioman_teleo.good.groups, dups=T)
summary.seqs(fasta=bioman_teleo.unique.good.filter.unique.pick.fasta, name=bioman_teleo.unique.good.filter.pick.names, processors=12)
count.groups(group=bioman_teleo.good.pick.groups)

#for clarity, rename files
system(cp bioman_teleo.unique.good.filter.unique.pick.fasta bioman_teleo_all.fasta)
system(cp bioman_teleo.unique.good.filter.pick.names bioman_teleo_all.names)
system(cp bioman_teleo.good.pick.groups bioman_teleo_all.groups)

#create count_table with the new files
make.table(name=bioman_teleo_all.names, group=bioman_teleo_all.groups, processors=12)

quit

# Clustering and taxonomic assignment

# VSEARCH + LULU
mothur "#cluster(fasta=bioman_teleo_all.fasta, count=bioman_teleo_all.count_table, method=agc)"
mv bioman_teleo_all.agc.unique_list.list bioman_teleo_all_vsearch.list
mothur "#get.oturep(list=bioman_teleo_all_vsearch.list, fasta=bioman_teleo_all.fasta, method=abundance, name=bioman_teleo_all.names)"
awk '/^>/{printf ">Otu%05d\n",++i; next}{print}' bioman_teleo_all_vsearch.0.03.rep.fasta | sed s/"\."//g| sed s/"-"//g > bioman_teleo_all_vsearch.0.03.rep_otuName.fasta 
mothur "#make.shared(list=bioman_teleo_all_vsearch.list, count=bioman_teleo_all.count_table, label=0.03)"
makeblastdb -in bioman_teleo_all_vsearch.0.03.rep_otuName.fasta -dbtype nucl
blastn -db bioman_teleo_all_vsearch.0.03.rep_otuName.fasta -outfmt '6 qseqid sseqid pident' -out bioman_teleo_all_vsearch.0.03.rep_otuName.dist -qcov_hsp_perc 100 -perc_identity 90 -query bioman_teleo_all_vsearch.0.03.rep_otuName.fasta -num_threads 8
cut -f 2,4- bioman_teleo_all_vsearch.shared > bioman_teleo_all_vsearch_0.03.tsv
Rscript LULU2.R bioman_teleo_all_vsearch.0.03.rep_otuName.dist bioman_teleo_all_vsearch_0.03.tsv 97
head -1 bioman_teleo_all_vsearch_0.03_curated.shared | cut -f 4- | sed s/"\t"/"\n"/g >  bioman_teleo_all_vsearch_0.03_curated.outs
nO=(`wc bioman_teleo_all_vsearch_0.03_curated.outs`)
n=`expr "$nO" + 3`
for f in `seq 4 $n`; do awk -v f="$f" '{s+=$f} END {print s}' bioman_teleo_all_vsearch_0.03_curated.shared; done > bioman_teleo_all_vsearch_0.03_curated.nRead_otu
grep -f bioman_teleo_all_vsearch_0.03_curated.outs -A1 bioman_teleo_all_vsearch.0.03.rep_otuName.fasta >  bioman_teleo_all_vsearch.0.03.rep_otuName_curated.fasta

mothur "#classify.seqs(fasta=bioman_teleo_all_vsearch.0.03.rep_otuName_curated.fasta, template=$DBlocation/$DBprefix.align, taxonomy=$DBlocation/$DBprefix.tax, method=wang, cutoff=30, processors=12)"
##mothur "#classify.seqs(fasta=bioman_teleo_all_vsearch.0.03.rep_otuName_curated.fasta, template=../database/teleo.align, taxonomy=../database/teleo.tax, method=wang, cutoff=30, processors=12)"

echo -e "OTU\tSize\tTaxonomy" > bioman_teleo_all_vsearch.0.03.rep_otuName_curated.rRNA_12S_teleo_worms_clean.wang_corr.taxonomy
paste <(cut -f1 bioman_teleo_all_vsearch.0.03.rep_otuName_curated.rRNA_12S_teleo_worms_clean.wang.taxonomy) <(cut -f1 bioman_teleo_all_vsearch_0.03_curated.nRead_otu) <(cut -f 2- bioman_teleo_all_vsearch.0.03.rep_otuName_curated.rRNA_12S_teleo_worms_clean.wang.taxonomy) | sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g >> bioman_teleo_all_vsearch.0.03.rep_otuName_curated.rRNA_12S_teleo_worms_clean.wang_corr.taxonomy


# SWARM + LULU
Rscript PREP_SWARM.R bioman_teleo_all.names bioman_teleo_all.fasta
swarm bioman_teleo_all_swarm.fasta -d 1 -t 10 -f -l bioman_teleo_all_swarm_d1.log -r -o bioman_teleo_all_swarm_d1.out -s bioman_teleo_all_swarm_d1.stats -w bioman_teleo_all_swarm_d1_rep.fasta 
awk '/^>/{printf ">Otu%04d\n",++i; next}{print}' bioman_teleo_all_swarm_d1_rep.fasta > bioman_teleo_all_swarm_d1_rep_otuName.fasta
sed -E  s/\(_[0-9]*,\)/","/g bioman_teleo_all_swarm_d1.out   | sed -E s/\(_[0-9]*\\t\)/"\\t"/g | sed -E  s/\(_[0-9]*$\)//g | sed s/swarm/swarm_1/g >  bioman_teleo_all_swarm_d1.list
mothur "#make.shared(list=bioman_teleo_all_swarm_d1.list, count=bioman_teleo_all.count_table, label=swarm_1)"
makeblastdb -in bioman_teleo_all_swarm_d1_rep_otuName.fasta -dbtype nucl
blastn -db bioman_teleo_all_swarm_d1_rep_otuName.fasta -outfmt '6 qseqid sseqid pident' -out bioman_teleo_all_swarm_d1_rep_otuName.dist -qcov_hsp_perc 100 -perc_identity 90 -query bioman_teleo_all_swarm_d1_rep_otuName.fasta -num_threads 8
cut -f 2,4- bioman_teleo_all_swarm_d1.shared > bioman_teleo_all_swarm_d1.tsv
Rscript LULU2.R bioman_teleo_all_swarm_d1_rep_otuName.dist bioman_teleo_all_swarm_d1.tsv 97
head -1 bioman_teleo_all_swarm_d1_curated.shared | cut -f 4- | sed s/"\t"/"\n"/g >  bioman_teleo_all_swarm_d1_curated.outs
nO=(`wc  bioman_teleo_all_swarm_d1_curated.outs`)
n=`expr "$nO" + 3`
for f in `seq 4 $n`; do awk -v f="$f" '{s+=$f} END {print s}'  bioman_teleo_all_swarm_d1_curated.shared; done >  bioman_teleo_all_swarm_d1_curated.nRead_otu
grep -f bioman_teleo_all_swarm_d1_curated.outs -A1 bioman_teleo_all_swarm_d1_rep_otuName.fasta >  bioman_teleo_swarm_d1_rep_otuName_curated.fasta

mothur "#classify.seqs(fasta=bioman_teleo_swarm_d1_rep_otuName_curated.fasta, template=$DBlocation/$DBprefix.align, taxonomy=$DBlocation/$DBprefix.tax, method=wang, cutoff=30, processors=12)"
echo -e "OTU\tSize\tTaxonomy" > bioman_teleo_swarm_d1_rep_otuName_curated.$DBprefix.wang_corr.taxonomy
paste <(cut -f1 bioman_teleo_swarm_d1_rep_otuName_curated.rRNA_12S_teleo_worms_clean.wang.taxonomy) <(cut -f1 bioman_teleo_all_swarm_d1_curated.nRead_otu) <(cut -f 2- bioman_teleo_swarm_d1_rep_otuName_curated.rRNA_12S_teleo_worms_clean.wang.taxonomy) | sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g >> bioman_teleo_swarm_d1_rep_otuName_curated.rRNA_12S_teleo_worms_clean.wang_corr.taxonomy
muscle
fasttree
	
# PHYLOTYPES (assign taxonomy without clustering)
mothur "#classify.seqs(fasta=bioman_teleo_all.fasta, template=$DBlocation/$DBprefix.align, taxonomy=$DBlocation/$DBprefix.tax, name=bioman_teleo_all.names, group=bioman_teleo_all.groups, method=wang, cutoff=30, processors=12)"
mothur "#phylotype(taxonomy=bioman_teleo_all.$DBprefix.wang.taxonomy)"
mothur "#make.shared(list=bioman_teleo_all.$DBprefix.wang.tx.list, count=bioman_teleo_all.count_table, label=1)"
mothur "#classify.otu(list=bioman_teleo_all.$DBprefix.wang.tx.list, count=bioman_teleo_all.count_table, taxonomy=bioman_teleo_all.$DBprefix.wang.taxonomy, label=1)"
sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g bioman_teleo_all.$DBprefix.wang.tx.1.cons.taxonomy > bioman_teleo_all.$DBprefix.wang.tx.1.cons_corr.taxonomy

# PHYLOTYPES GLOBAL DB (assign taxonomy without clustering)
mothur "#classify.seqs(fasta=bioman_teleo_all.fasta, template=$globalDBlocation/$globalDBprefix.align, taxonomy=$globalDBlocation/$globalDBprefix.tax, name=bioman_teleo_all.names, group=bioman_teleo_all.groups, method=wang, cutoff=30, processors=12)"
mothur "#phylotype(taxonomy=bioman_teleo_all.$globalDBprefix.wang.taxonomy)"
mothur "#make.shared(list=bioman_teleo_all.$globalDBprefix.wang.tx.list, count=bioman_teleo_all.count_table, label=1)"
mothur "#classify.otu(list=bioman_teleo_all.$globalDBprefix.wang.tx.list, count=bioman_teleo_all.count_table, taxonomy=bioman_teleo_all.$globalDBprefix.wang.taxonomy, label=1)"
sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g bioman_teleo_all.$globalDBprefix.wang.tx.1.cons.taxonomy > bioman_teleo_all.$globalDBprefix.wang.tx.1.cons_corr.taxonomy