/*  array.h -- simple array handling functions.

    Copyright (C) 2016 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __ARRAY_H__
#define __ARRAY_H__

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
 * integer array functions
 */

typedef struct {
    int end;
    int max;
    int *entries;
} ia_t;

ia_t *ia_init(int max);
int ia_compare(const void *ia1, const void *ia2);
void ia_sort(ia_t *ia);
void ia_push(ia_t *ia, int i);
void ia_free(ia_t *ia);
char *ia_join(ia_t *ia, char *delim);
int ia_sum(ia_t *ia);
static inline bool ia_isEmpty(ia_t *ia) { return (ia->end == 0); }


/*
 * generic arrays
 */

typedef struct {
    int end;
    int max;
    void (*free_entry)(void *);
    void **entries;
} va_t;

va_t *va_init(int max, void(*free_entry)(void*));
void va_push(va_t *va, void *ent);
void va_free(va_t *va);
static inline bool va_isEmpty(va_t *va) { return va->end == 0; }

// NB these functions assume an array of strings
char *va_join(va_t *va, char *delim);
int va_contains(va_t *va, char *s);

#endif

