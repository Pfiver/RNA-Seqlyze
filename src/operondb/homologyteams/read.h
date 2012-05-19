/****************************************************************
 * This file has been modified from its original version -
 * homologyteams, version 1.0, provided by Michael H. Goldwasser 
 * (goldwamh@slu.edu) and Xin He (xinhe2@uiuc.edu)
 * 
 * Copyright (C) 2008  Mihaela Pertea (mpertea@umiacs.umd.edu)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *****************************************************************
 *
 * File:    read.h
 *
 * Support for reading files.
 *   -- file format for original definition of chromosomes pair
 *   -- file format for delimited team reports
 *
 ****************************************************************/
#ifndef READ_H
#define READ_H

#include <stdio.h>
#include "types.h"


#define MAX_GENE_NAME 20
#define MAX_LINE 100000
#define MAX_DIGITS 20

/*
 * read the original chromosome file, constructing both chromosomes, as well as the
 * cog_names table.
 *
 * returns 0 if success, -1 if failure
 */
int readChromosomes(FILE* data_file, Chromosome* chrC, Chromosome* chrD, char*** family_names );


/*
 * Reads next delimiter from the given file, and sets the given delim
 * appropriately.
 *
 * Returns true if success, false if failure
 */
boole readNextDelim(FILE* delimFile, Delimiter *delim);


#endif

