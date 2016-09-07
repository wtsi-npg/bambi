/*  test/i2b/filterfile.c -- filterfile test cases.

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

#include "../../src/filterfile.c"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

int verbose = 0;

int success = 0;
int failure = 0;

void checkLike(char *name, char *expected, char *actual)
{
    if (actual == NULL) actual = "<null>";
    if (strstr(actual, expected) == NULL) {
        fprintf(stderr, "%s\n" "Expected: %s\n" "Got:      %s\n", name, expected, actual);
        failure++;
    }
}

void checkEqual(char *name, char *expected, char *actual)
{
    if (actual == NULL) actual = "<null>";
    if (strcmp(expected, actual)) {
        fprintf(stderr, "%s\n" "Expected: %s\n" "Got:      %s\n", name, expected, actual);
        failure++;
    }
}

void icheckEqual(char *name, int expected, int actual)
{
    if (expected != actual) {
        fprintf(stderr, "%s\n" "Expected: %d\n" "Got:      %d\n", name, expected, actual);
        failure++;
    }
}

int main(int argc, char**argv)
{
    int n, f;

    filter_t *filter;

    filter = filter_open(MKNAME(DATA_DIR,"/s_1_1101.filter"));
    if (filter->errmsg) {
        fprintf(stderr,"Error opening file '%s':  %s\n", MKNAME(DATA_DIR,"/s_1_1101.filter"), filter->errmsg);
        failure++;
    }
    icheckEqual("Version", 3, filter->version);
    icheckEqual("Total clusters", 2609912, filter->total_clusters);
    icheckEqual("Current cluster", 0, filter->current_cluster);

    f = filter_next(filter);
    icheckEqual("Next Current cluster", 1, filter->current_cluster);
    icheckEqual("Next Current PF clusters", 0, filter->current_pf_cluster);
    icheckEqual("First entry", 0, f);

    for (n=0; n<318; n++) {
        f = filter_next(filter);
    }
    icheckEqual("319 entry", 1, f);
    icheckEqual("319 Current cluster", 319, filter->current_cluster);
    icheckEqual("319 Current PF clusters", 264, filter->current_pf_cluster);
    icheckEqual("319 Total clusters", 2609912, filter->total_clusters);

    while (filter_next(filter) != -1);
    icheckEqual("Last Current cluster", 2609912, filter->current_cluster);
    icheckEqual("Last Current PF clusters", 2425954, filter->current_pf_cluster);
    icheckEqual("Last Total clusters", 2609912, filter->total_clusters);

    printf("filter tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
