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
 * File:    filter.h
 *
 * This file includes a variety of filters which control the precise
 * rules used in selecting which teams are included in the output.
 *
 ****************************************************************/
#ifndef FILTER_H
#define FILTER_H

#include "globaldata.h"
#include "types.h"


typedef struct filter_ {

  /**
   * identifier used in matching command-line argument
   */
  char* name;                                

  /**
   * check parameter
   */
  boole (*checkParam)();

  /**
   * called once, during initialization of algorithm
   */
  void (*begin)();

  /**
   * called at the time each team is first identified
   */
  void (*reportTeam)(Delimiter);

  /**
   * called once, after the algorithm completes
   */
  void (*end)();

  /**
   * description of the filter
   */
  char* description;                                

  /**
   * description of the parameter
   */
  char* paramdesc;

} filter;



/*
 *
 * The array of possible filters from which to choose
 *
 */
extern filter filters[];

/*
 *
 * The index of the chosen filter
 *
 */
extern int chosenFilter;


/*
 *  Sets the chosenFilter variable, based on the given 'key'.  For some
 *  of the choices, param allows for specification of an optional
 *  parameter.
 *
 *  returns failure if key not found
 */
boole chooseFilter(char* key, char* param);


/*
 *  Prints the keywords and descriptions for the existing filters
 */
void printFilterChoices(char* preface);


#endif
