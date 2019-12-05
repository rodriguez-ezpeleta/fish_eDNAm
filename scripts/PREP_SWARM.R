setwd(".")

args <- commandArgs()
names <- args[6]
fasta <- args[7]
swarm_fasta <- gsub(".fasta", "_swarm.fasta", fasta)

names_file <- scan(file=names, what="", quiet=TRUE) 
names_data <- names_file[(1:length(names_file)) %% 2 == 0] 
names(names_data) <- names_file[(1:length(names_file)) %% 2 == 1] 
n_seqs <- nchar(names_data) - nchar(gsub(",", "", names_data)) + 1 
  
  
fasta_data <- scan(fasta, what="", quiet=TRUE) 
sequence_data <- fasta_data[grepl("^[ATGCatgc.-]", fasta_data)] 
sequence_data <- gsub("[-.]", "", sequence_data) 
names(sequence_data) <- gsub(">", "", fasta_data[grepl("^>", fasta_data)], 2, )   
  
seq_with_freq <- paste0(">", names(sequence_data), "_", n_seqs[names(sequence_data)], "\n", sequence_data) 
    
write(seq_with_freq, swarm_fasta) 
  
