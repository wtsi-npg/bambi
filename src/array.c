/*  array.c -- simple array handling functions.

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

#include "array.h"
#include <assert.h>

/*
 * integer array functions
 */

int ia_compare(const void *ia1, const void *ia2)
{
    return *(int *)ia1 - *(int *)ia2;
}

void ia_sort(ia_t *ia)
{
    qsort(ia->entries, ia->end, sizeof(int), ia_compare);
}

void ia_push(ia_t *ia, int i)
{
    if (ia->end == ia->max) {
        // expand the array
        ia->max *= 2;
        ia->entries = realloc(ia->entries, ia->max * sizeof(int));
    }
    ia->entries[ia->end] = i;
    ia->end++;
}

void ia_free(ia_t *ia)
{
    if (!ia) return;
    free(ia->entries);
    free(ia);
}

ia_t *ia_init(int max)
{
    ia_t *ia = calloc(1, sizeof(ia_t));
    ia->end = 0;
    ia->max = max;
    ia->entries = calloc(ia->max, sizeof(int));
    return ia;
}

char *ia_join(ia_t *ia, char *delim)
{
    int m = 64;
    char *s = calloc(m,1);
    char *a = calloc(64,1);

    for (int n=0; n < ia->end; n++) {
        while (strlen(s)+strlen(delim)>=m) { m *= 2; s = realloc(s,m); }
        if (n) strcat(s,delim);
        sprintf(a,"%d",ia->entries[n]);
        assert(strlen(a)<64);
        while (strlen(s)+strlen(a)>=m) { m *= 2; s = realloc(s,m); }
        strcat(s,a);
        assert(strlen(s) < m);
    }
    free(a);
    return s;
}

int ia_sum(ia_t *ia)
{
    int sum=0;
    for (int n=0; n < ia->end; n++) {
        sum += ia->entries[n];
    }
    return sum;
}


/*
 * generic arrays
 */

va_t *va_init(int max, void(*free_entry)(void*))
{
    va_t *va = calloc(1,sizeof(va_t));
    va->end = 0;
    va->max = max;
    va->free_entry = free_entry;
    va->entries = calloc(va->max, sizeof(void *));
    return va;
}

void va_push(va_t *va, void *ent)
{
    if (va->end == va->max) {
        // expand the array
        va->max *= 2;
        va->entries = realloc(va->entries, va->max * sizeof(void *));
    }
    va->entries[va->end] = ent;
    va->end++;
}

void va_free(va_t *va)
{
    if (!va) return;
    if (va->free_entry != NULL) {
        for (int n=0; n < va->end; n++) {
            va->free_entry(va->entries[n]);
        }
    }
    free(va->entries);
    free(va);
}

/*
 * join entries into a 'delim' separated list
 * NB assumes that entries are strings
 */
char *va_join(va_t *va, char *delim)
{
    int m = 64;
    char *s = calloc(m,1);
    char *a;

    for (int n=0; n < va->end; n++) {
        while (strlen(s)+strlen(delim)>=m) { m *= 2; s = realloc(s,m); }
        if (n) strcat(s,delim);
        a = va->entries[n];
        while (strlen(s)+strlen(a)>=m) { m *= 2; s = realloc(s,m); }
        strcat(s,a);
        assert(strlen(s) < m);
    }
    return s;
}

/*
 * return array number containing 's', or -1 if not found
 * NB assumes all entries are strings
 */
int va_contains(va_t *va, char *s)
{
    for (int n=0; n < va->end; n++) {
        if (strcmp(va->entries[n],s) == 0) return n;
    }
    return -1;
}

