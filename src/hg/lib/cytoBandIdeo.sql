# cytoBand.sql was originally generated by the autoSql program, which also 
# generated cytoBand.c and cytoBand.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#Describes the positions of cytogenetic bands with a chromosome
CREATE TABLE cytoBandIdeo (
    chrom varchar(255) not null,	# Human chromosome number
    chromStart int unsigned not null,	# Start position in genoSeq
    chromEnd int unsigned not null,	# End position in genoSeq
    name varchar(255) not null,	# Name of cytogenetic band
    gieStain varchar(255) not null,	# Giemsa stain results
              #Indices
    PRIMARY KEY(chrom(12),chromStart),
    UNIQUE(chrom(12),chromEnd)
);
