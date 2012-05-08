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
 * File: types.h
 *
 * This file includes the type definitions for genes, chromosomes,
 * subchromosome, and team delimeters.
 *
 ****************************************************************/

#ifndef TYPES_H
#define TYPES_H

/*
 * representation of a single gene, occuring within a chromosome
 */
typedef enum {F=0, T=1} boole;

/*
 * representation of a single gene, occuring within a chromosome
 */
typedef struct Gene_ {
  int	 family;		// The containing family for this gene (based on unique FAMILY ID)
  // before: float position;      // The gene's position within its chromosome
  long positionBeg;      // The gene's begin position within its chromosome; [MP] code
  long positionEnd;      // The gene's end position within its chromosome; [MP] code
  int geneord;         // The gene's order within its chromosome; [MP] code
} Gene;


/*
 * representation of a chromosome
 */
typedef struct Chromosome_ {
  char  *name;	    // the name of tag of a chromosome
  int   numGenes;   // number of genes represented in this chromosome
  Gene  *genes;     // the array of genes
} Chromosome;


/*
 * compact representation of a subchromosome
 */
typedef struct Subchrom_ {
  Chromosome*	chrom;    
  int		startIndex;	// note that demarcation is by index, not position
  int		endIndex;       // note that demarcation is by index, not position
} Subchrom;


/* 
 * compact representation of a COG team 
 */
typedef struct Delimiter_ {
  Subchrom  subC;
  Subchrom  subD;
} Delimiter; 


#endif
