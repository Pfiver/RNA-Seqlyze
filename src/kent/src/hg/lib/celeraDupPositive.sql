# first.sql was originally generated by the autoSql program, which also 
# generated first.c and first.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#Summary of large genomic Duplications from Celera Data
CREATE TABLE celeraDupPositive (
    bin smallint not null,      # Index field
    chrom varchar(255) not null,	# Human chromosome or FPC contig
    chromStart int unsigned not null,	# Start position in chromosome
    chromEnd int unsigned not null,	# End position in chromosome
    name varchar(255) not null,	# Celera accession Name
    fullname varchar(255) not null,	# Celera Accession Fullname
    fracMatch float not null,	# fraction of matching bases
    bpalign float not null,	# base pair alignment score
              #Indices
    INDEX (chrom(8), bin),
    INDEX (chrom(8), chromStart),
    INDEX (name(18))
);
