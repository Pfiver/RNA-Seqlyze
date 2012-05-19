# allenBrainUrl.sql was originally generated by the autoSql program, which also 
# generated allenBrainUrl.c and allenBrainUrl.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#Link together allenBrain ID and a URL
CREATE TABLE allenBrainUrl (
    name varchar(255) not null,	# Allen Brain Atlas ID
    url varchar(255) not null,	# URL of link into Allen site
              #Indices
    PRIMARY KEY(name)
);
