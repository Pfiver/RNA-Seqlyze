# rnaGroup.sql was originally generated by the autoSql program, which also 
# generated rnaGroup.c and rnaGroup.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#Details of grouping of rna in cluster
CREATE TABLE rnaGroup (
    chrom varchar(255) not null,	# Chomosome
    chromStart int unsigned not null,	# Start position in chromosome
    chromEnd int unsigned not null,	# End position in chromosome
    name varchar(255) not null,	# Name of item
    score int unsigned not null,	# Always 1000
    strand char(1) not null,	# + or -
    refSeqCount int unsigned not null,	# Number of refSeq mrnas in cluster
    refSeqs longblob not null,	# List of refSeq accessions
    genBankCount int unsigned not null,	# Number of Genbank mrnas in cluster
    genBanks longblob not null,	# List of Genbank accessions
    rikenCount int unsigned not null,	# Number of Riken mrnas in cluster
    rikens longblob not null,	# List of Riken IDs
              #Indices
    PRIMARY KEY(name)
);
