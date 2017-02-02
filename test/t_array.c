/*  t_array.c -- array test cases.

    Copyright (C) 2016 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "array.h"
#include "htslib/hts.h"

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)


int verbose = 0;

int success = 0;
int failure = 0;

void checkEqual(char *name, char *expected, char *actual)
{
    if (actual == NULL) actual = "<null>";
    if (strcmp(expected, actual)) {
        fprintf(stderr, "%s: Expected: %s \tGot: %s\n", name, expected, actual);
        failure++;
    }
}

void icheckEqual(char *name, int expected, int actual)
{
    if (expected != actual) {
        fprintf(stderr, "%s: Expected: %d \tGot: %d\n", name, expected, actual);
        failure++;
    }
}

int main(int argc, char**argv)
{
    va_t *va = va_init(2,free);
    ia_t *ia = ia_init(3);

    icheckEqual("va Initially empty", true, va_isEmpty(va));
    icheckEqual("ia Initially empty", true, ia_isEmpty(ia));

    va_push(va,strdup("Hello world"));
    ia_push(ia, 1);
    icheckEqual("va one item", 1, va->end);    
    icheckEqual("ia one item", 1, ia->end);    

    va_push(va,strdup("xyzzy"));
    ia_push(ia, 2);
    va_push(va,strdup("xyzzy"));
    ia_push(ia, 3);
    va_push(va,strdup("xyzzy"));
    ia_push(ia, 9);
    va_push(va,strdup("plugh"));
    ia_push(ia, 5);
    icheckEqual("va five items", 5, va->end);    
    icheckEqual("ia five items", 5, ia->end);    

    checkEqual("va first item", "Hello world", va->entries[0]);
    icheckEqual("ia first item", 1, ia->entries[0]);
    checkEqual("va last item", "plugh", va->entries[va->end-1]);
    icheckEqual("ia last item", 5, ia->entries[ia->end-1]);

    ia_sort(ia);
    icheckEqual("ia sort first", 1, ia->entries[0]);
    icheckEqual("ia sort last", 9, ia->entries[ia->end-1]);

    char *s = ia_join(ia,"xyz");
    checkEqual("ia join", "1xyz2xyz3xyz5xyz9", s);

    icheckEqual("va_contains(1)", 1, va_contains(va,"xyzzy"));
    icheckEqual("va_contains(0)", 0, va_contains(va,"Hello world"));
    icheckEqual("va_contains(4)", 4, va_contains(va,"plugh"));
    icheckEqual("va_contains(-1)", -1, va_contains(va,"Garp"));

    ia_free(ia);
    va_free(va);

    printf("array tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
