# wabaCrude.sql was originally generated by the autoSql program, which also 
# generated wabaCrude.c and wabaCrude.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#Info on a crude alignment
CREATE TABLE wabaCrude (
    score int unsigned not null,	# score 0 to about 6000
    qFile varchar(255) not null,	# query sequence file
    qSeq varchar(255) not null,	# query sequence within file
    qStart int unsigned not null,	# query start position
    qEnd int unsigned not null,	# query end position
    strand tinyint not null,	# strand of alignment: +1 or -1
    tFile varchar(255) not null,	# target sequence file
    tSeq varchar(255) not null,	# target sequence within file
    tStart int unsigned not null,	# target start position
    tEnd int unsigned not null,	# target end position
              #Indices
    PRIMARY KEY(score)
);
