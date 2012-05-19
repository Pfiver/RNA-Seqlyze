# ntOoaHaplo.sql was originally generated by the autoSql program, which also 
# generated ntOoaHaplo.c and ntOoaHaplo.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#Candidate regions for gene flow from Neandertal to non-African modern humans (Table 5 of Green RE et al., Science 2010)
CREATE TABLE ntOoaHaplo (
    bin smallint(5) unsigned not null,  # index for browser speedup
    chrom varchar(255) not null,	# Reference sequence chromosome
    chromStart int unsigned not null,	# Start position in chromosome
    chromEnd int unsigned not null,	# End position in chromosome
    name varchar(255) not null,	# Qualitative assessment (OOA = out of africa, COS = cosmopolitan)
    score int unsigned not null,	# For BED compatibility: Score from 0-1000 (placeholder = 0)
    strand char(1) not null,	# For BED compatibility: + or - (placeholder = +)
    thickStart int unsigned not null,	# For BED compatibility: Start of where display should be thick
    thickEnd int unsigned not null,	# For BED compatibility: End of where display should be thick
    reserved int unsigned not null,	# itemRgb color code
    st float not null,	# Estimated ratio of OOA/African gene tree depth
    ooaTagFreq float not null,	# Average frequency of tag in OOA clade
    am tinyint unsigned not null,	# Neandertal (M)atches OOA-specific clade (Ancestral)
    dm tinyint unsigned not null,	# Neandertal (M)atches OOA-specific clade (Derived)
    an tinyint unsigned not null,	# Neandertal does (N)ot match OOA-specific clade (Ancestral)
    dn tinyint unsigned not null,	# Neandertal does (N)ot match OOA-specific clade (Derived)
              #Indices
    INDEX chrom (chrom,bin)
);