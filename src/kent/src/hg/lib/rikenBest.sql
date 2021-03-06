# rikenBest.sql was originally generated by the autoSql program, which also 
# generated rikenBest.c and rikenBest.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#The best Riken mRNA in each cluster
CREATE TABLE rikenBest (
    name varchar(255) not null,	# Riken mRNA ID
    orfScore float not null,	# Score from bestorf - 50 is decent
    orfStrand char(1) not null,	# Strand bestorf is on: +, - or .
    intronOrientation int not null,	# +1 for each GT/AG intron, -1 for each CT/AC
    position varchar(255) not null,	# Position in genome of cluster chrN:start-end format
    rikenCount int not null,	# Number of Riken mRNAs in cluster
    genBankCount int not null,	# Number of Genbank mRNAs in cluster
    refSeqCount int not null,	# Number of RefSeq mRNAs in cluster
    clusterId varchar(255) not null,	# ID of cluster
              #Indices
    PRIMARY KEY(clusterId(16)),
    INDEX(name(12))
);
