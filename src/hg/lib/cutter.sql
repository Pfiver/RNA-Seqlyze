# cutter.sql was originally generated by the autoSql program, which also 
# generated cutter.c and cutter.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#Restriction Enzymes
CREATE TABLE cutter (
    name varchar(255) not null,	# Name of Enzyme
    size int unsigned not null,	# Size of recognition sequence
    matchSize int unsigned not null,	# size without N's
    seq varchar(255) not null,	# Recognition sequence
    cut int unsigned not null,	# Cut site on the plus strand
    overhang int not null,	# Overhang size and direction
    palindromic tinyint unsigned not null,	# 1 if it's panlidromic, 0 if not.
    semicolon tinyint unsigned not null,	# REBASE record: 0 if primary isomer, 1 if not.
    numSciz int unsigned not null,	# Number of isoschizomers
    scizs longblob not null,	# Names of isoschizomers
    numCompanies int unsigned not null,	# Number of companies selling this enzyme
    companies longblob not null,	# Company letters
    numRefs int unsigned not null,	# Number of references
    refs longblob not null	# Reference numbers
              #Indices
);
