/****************************************************************
 *
 * Project: homologyteams
 * Authors: Michael H. Goldwasser (goldwamh@slu.edu) and Xin He (xinhe2@uiuc.edu)
 * Version: 1.0 (May 2004)
 *
 * Copyright (C) 2004  Michael H. Goldwasser and Xin He
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
 * File:    output.h
 *
 * Support for outputting a desired team, in an appropriate format,
 * where a team is identified in delimited form based upon the two
 * identifying subchromosomes.  Currently, there exists support for
 * three formats:
 *     DELIM  - abbreviated output of original delimters
 *     COGS   - output the COGs
 *     WITNESS - output the COGs, as well as the witnesses
 *
 ****************************************************************/
#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>
#include "globaldata.h"
#include "types.h"


/*
 *  The output file.  Will be stdout by default, but can be reset
 *  based upon command-line arguments.
 */ 
extern FILE* ofp;

/*
 *  Generates output for the team defined by the give delimiter.  
 *  The style of the output is based upon the current setting of
 *  outputStyle.
 *
 *  If not NULL, the second parameter is an array of the COG IDs,
 *  with the third parameter the size of the array.
 */ 
void outputTeam(Delimiter delim, int** knownContents, int knownSize);


/*
 *  Sets the desired output style, based on the given 'key'.
 *
 *  returns failure if key not found
 */
boole chooseOutputStyle(char* key);


/*
 *  Prints the keywords for the existing output choices
 */
void printOutputStyles(char* preface);



#endif
