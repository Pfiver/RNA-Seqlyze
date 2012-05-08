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
 * File:    list.h
 *
 * This file includes the definitions for a singly-linked list, using
 * a header sentinal.
 *
 ****************************************************************/

#ifndef LIST_H
#define LIST_H

typedef struct Node_ {
  struct Node_* next;
  void* item;
} Node;

typedef struct List_ {
  Node* head;
  Node* last;
  int   size;
} List;



/*
 *  Creates and returns a new (empty) list
 */
List* listNew();

/*
 *  Return size of given list (not including the sentinal)
 */
int listSize(List* thisL);

/*
 *  Deallocate the given list
 *  (though not the elements contained)
 */
void listFree(List* thisL);

/*
 *  Add the new element to the front of the list
 */
void listAddFront(List* thisL, void* element);

/*
 *  Add the new element to the back of the list
 */
void listAddBack(List* thisL, void* element);



#endif
