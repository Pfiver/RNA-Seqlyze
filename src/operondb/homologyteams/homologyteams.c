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
 * File:    homologyteams.c
 *
 * This provides the main driver for the software.
 *
 ****************************************************************/
#include <stdlib.h>
#include <string.h>
#include "algorithm.h"
#include "filter.h"
#include "globaldata.h"
#include "output.h"
#include "read.h"


/****************************************************************
 *
 * True definition of the global data declared in
 *    globaldata.h
 *
 ****************************************************************/

// before: double delta = 1000.0;  // the parameter for maximum allowable gap in position
long delta = 1000.0;  // the parameter for maximum allowable gap in position; [MP] : changed to long since I don't see why it should be double
int no_genes=1;         // [MP] distance measured in no of genes 
int   num_families;	        // the total number of homology families
char  **family_names;      // array of the family names
int bparam=0;           // [MP] how to compute the distance between genes:
                        // b=0 : use only delta parameter (default)
                        // b=1 : use the no_genes parameters
                        // b=2 : use both parameters
int circC=0;         // [MP] 0 if chr C linear; 1 if circular (should be the default)
int circD=0;         // [MP] 0 if chr D linear; 1 if circular (should be the default)

/****************************************************************
 *
 * Variables local to this file
 *
 ****************************************************************/

FILE  *inputfp = NULL;
FILE  *prefp = NULL;
Chromosome C;
Chromosome D;

/****************************************************************
 *
 * Some utility methods
 *
 ****************************************************************/


void printUsage () 
{
  printf("\n");
  printf(" Usage:\n");
  printf("   homologyteams [options] inputfile\n");
  printf("\n");
  printf(" Options:\n");
  printf("   -h       this message\n");
  printf("\n");
  printf("   -d val   delta value (default is 1000.0)\n");
  printf("\n");
  printf("   -n val   val gives distance in number of genes\n"); // start [MP] code
  printf("\n");
  printf("   -b val   val=0 : use delta parameter to compute distance (default);\n");
  printf("            val=1 : use number of genes to compute distance;\n");
  printf("            val=2 : use both measures to compute distance;\n");
  printf("\n"); // end [MP] code
  printf("   -W file  write output to the given file\n");
  printf("\n");
  printf("   -F rule  use filtering rule for selecting relevant teams\n");
  printf("   -P param additional parameter for chosen filtering rule\n");
  printf("\n");
  printFilterChoices("            ");
  printf("\n");
  printf("   -O style use given output style\n");
  printOutputStyles("            ");
  printf("\n");
  printf("   -T file  use teams from precomputed file rather than\n");
  printf("              computing teams from scratch (see -O concise option)\n");
  printf("\n");
  /* Just don't use this for now
  printf("   -c val   val=0 if both chromosomes are linear (default);\n"); // start [MP] code
  printf("            val=1 if first chromosome only is circular");
  printf("            val=2 if second chromosome only is circular");
  printf("            val=3 if both chromosomes are circular"); //should be made default
  printf("\n"); // end [MP] code*/

}



int processArguments (int argc, char* argv[]) {
  int failure = F;
  int i;
  int circular; 

  char* filterName = NULL;
  char* filterParam = NULL;
  ofp = stdout;                   // Default

  for (i=1; !failure && i<argc; i++) {
    if (strcmp(argv[i],"-h")==0) {
      failure = T;
    } else if (strcmp(argv[i],"-d")==0) {
      if (argc<i+1) {
	fprintf(stderr,"delta value must be specified after -d option\n");
	failure = T;
      } else {
	char* temp;
	delta = strtol(argv[++i],&temp,10);
	if (delta <= 0 || temp==argv[i]) {
	  // ERROR occurred
	  fprintf(stderr,"Specified delta value (%s) is invalid.\n",argv[i]);
	  failure = T;
	}
      }
    } else if (strcmp(argv[i],"-c")==0) {  // start [MP] code here
      if(argc<i+1) {
	fprintf(stderr,"circularity value must be specified after -c option\n");
	failure = T;
      } else {
	circular = atoi(argv[++i]);
	if(circular !=0 && circular !=1 && circular !=2 && circular !=3) {
	  fprintf(stderr,"Circular value must be {0,1,2,3}.\n");
	  failure = T;
	}
	switch (circular) {
	case 0: circC=0; circD=0;break;
	case 1: circC=1; circD=0;break;
	case 2: circC=0; circD=1;break;
	case 3: circC=0; circD=1;break;
	}
      }
    } else if(strcmp(argv[i],"-n")==0) { 
      if(argc<i+1) {
	fprintf(stderr,"Number of genes must be specified after -n option\n");
	failure = T;
      } else {
	no_genes = atoi(argv[++i]);
	if(no_genes<1) {
	  fprintf(stderr,"Number of genes must be strict positive.\n");
	  failure = T;
	}
      }
    } else if(strcmp(argv[i],"-b")==0) { 
      if(argc<i+1) {
	fprintf(stderr,"Incorrect value after -b option\n");
	failure = T;
      } else {
	bparam = atoi(argv[++i]);
	if(bparam !=0 && bparam != 1 && bparam !=2) {
	  fprintf(stderr,"Specified value %d after -b option is invalid.\n",bparam);
	  failure = T;
	}
      } // end [MP] code
    } else if (strcmp(argv[i],"-W")==0) {
      if (argc<i+1) {
	fprintf(stderr,"desired filename must be specified after -W option\n");
	failure = T;
      } else if ((ofp=fopen(argv[++i],"w"))==NULL) {
	fprintf(stderr,"Could not open output file %s.\n",argv[i]);
	failure = T;
      }
    } else if (strcmp(argv[i],"-T")==0) {
      if (argc<i+1) {
	fprintf(stderr,"precomputed file of teams must be specified after -T option\n");
	failure = T;
      } else if ((prefp=fopen(argv[++i],"r"))==NULL) {
	fprintf(stderr,"Could not read precomuted file of teams (%s).\n",argv[i]);
	failure = T;
      }
    } else if (strcmp(argv[i],"-O")==0) {
      if (argc<i+1) {
	fprintf(stderr,"output style must be specified after -O option\n");
	failure = T;
      } else {
	if (!chooseOutputStyle(argv[++i]))
	  failure = T;
      }
    } else if (strcmp(argv[i],"-F")==0) {
      if (argc<i+1) {
	fprintf(stderr,"filter rule must be specified after -F option\n");
	failure = T;
      } else {
	filterName = argv[++i];
      }
    } else if (strcmp(argv[i],"-P")==0) {
      if (argc<i+1) {
	fprintf(stderr,"parameter string must be specified after -P option\n");
	failure = T;
      } else {
	filterParam = argv[++i];
      }
    } else if (inputfp==NULL) {
      if ((inputfp=fopen(argv[i],"r"))==NULL) {
	fprintf(stderr,"Could not read input file %s.\n",argv[i]);
	failure = T;
      }
    } else {
      fprintf(stderr,"Unrecognized argument (%s)\n",argv[i]);
      failure = T;
    }
  }

  if (!failure) {
    if (!inputfp) {  // required argument not given
      fprintf(stderr,"You must specify and input file contaiing two chromosomes.\n");
      failure = T;
    } else if (!prefp && (delta == -1)) {     // required arguments not reached
      fprintf(stderr,"delta is a required argument.\n");
      failure = T;
    } else if (filterName && !chooseFilter(filterName,filterParam)) {
      failure=T;
    }
  }

  if (prefp && (delta!=1000.0))
    fprintf(stderr,"WARNING:  specified delta value is ignored when using precomputed teams.\n");


  if (failure)
    printUsage();

  return failure;
}




int main (int argc, char* argv[]) 
{     
  int failure;


  if ((failure=processArguments(argc,argv))==0) {

    if (readChromosomes(inputfp,&C,&D,&family_names)!=0) {
      fprintf(stderr,"Unable to properly read the original input file.  Please check the format of the data.\n");
      failure=1;
    }
    fclose(inputfp);



    if (!failure) {
      if (prefp) {
	// use precalculated delimeter file to identify teams
	Delimiter delim;
	delim.subC.chrom = &C;
	delim.subD.chrom = &D;

	// since teams are begin processed without algorithm, must
        // expliclity initialize the necessary structures for the
        // algorithm module.
	initializeTime();

	filters[chosenFilter].begin();
	while (readNextDelim(prefp,&delim))
	  filters[chosenFilter].reportTeam(delim);
	filters[chosenFilter].end();

	// and explicitly free the data created by initializetime()
	freeGlobal();


      } else {
	// run the algorithm to identify teams
	findTeams(&C,&D);
    
      }

    }


    // deallocate various data
    if (C.genes) free(C.genes);
    if (D.genes) free(D.genes);
    if (family_names) {
      int i;
      for (i=0; i<num_families; i++)
	if (family_names[i]) free(family_names[i]);
    }
    free(family_names);
  }

  if ((ofp != NULL) && (ofp !=stdout))
    fclose(ofp);

  return failure;
}

