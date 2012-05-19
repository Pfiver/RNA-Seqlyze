# dgv.sql was originally generated by the autoSql program, which also 
# generated dgv.c and dgv.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#Database of Genomic Variants
CREATE TABLE dgv (
    bin smallint not null,		# Bin number for browser speedup
    chrom varchar(255) not null,	# Reference sequence chromosome or scaffold
    chromStart int unsigned not null,	# Start position in chromosome
    chromEnd int unsigned not null,	# End position in chromosome
    name varchar(255) not null,	# Name of item
    score int unsigned not null,	# Score from 0-1000
    strand char(1) not null,	# + or -
    thickStart int unsigned not null,	# Start of where display should be thick (start codon)
    thickEnd int unsigned not null,	# End of where display should be thick (stop codon)
    itemRgb int unsigned not null,	# Item R,G,B color.
    landmark varchar(255) not null,	# Genomic marker near the variation locus
    varType varchar(255) not null,	# Type of variation
    reference varchar(255) not null,	# Literature reference for the study that included this variant
    pubMedId int unsigned not null,	# For linking to pubMed abstract of reference
    method varchar(255) not null,	# Brief description of method/platform
    sample longblob not null,	# Description of sample population for the study
              #Indices
    INDEX (chrom,bin),
    PRIMARY KEY (name)
);