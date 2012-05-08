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
 * File:    globaldata.h
 *
 * Declares several pieces of data which are defined in homologyteams.c but
 * which need to be accessed from within several other modules.
 *
 ****************************************************************/
#ifndef GLOBALDATA_H
#define GLOBALDATA_H

// before: extern double delta;    // the parameter for maximum allowable gap in position
extern long delta;    // the parameter for maximum allowable gap in position; [MP] : changed to long
extern int no_genes;    // [MP] distance measured in no of genes 
extern bparam;          // [MP] how to compute the distance between genes:
                        // b=0 : use only delta parameter
                        // b=1 : use the no_genes parameters
                        // b=2 : use both parameters
extern int circC;         // [MP] c=1 if chromosome is circular (should be made default)
extern int circD;         // [MP] c=1 if chromosome is circular (should be made default)
extern int   num_families;    // the total number of COGs
extern char **family_names;  // array of the COG names

#endif
