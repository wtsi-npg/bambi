/*  test/i2b/posfile.c -- posfile test cases.

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

#include "posfile.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

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
    int result, n;

    posfile_t *posfile;

    posfile = posfile_open("test/i2b/110323_HS13_06000_B_B039WABXX/Data/Intensities/L001/s_1_1101.clocs");
    if (posfile->errmsg) {
        fprintf(stderr,"Error opening file: %s\n", posfile->errmsg);
        failure++;
    }
    icheckEqual("Version", 1, posfile->version);
    icheckEqual("Total blocks", 65600, posfile->total_blocks);

    posfile_next(posfile);
    icheckEqual("next X", 1235, posfile_get_x(posfile));
    icheckEqual("next Y", 1989, posfile_get_y(posfile));
    icheckEqual("current block", 247, posfile->current_block);

    for (n=0; n<306; n++) {
        posfile_next(posfile);
    }
    icheckEqual("307 x", 1279, posfile_get_x(posfile));
    icheckEqual("307 y", 2120, posfile_get_y(posfile));
    icheckEqual("307  block", 330, posfile->current_block);

    while (posfile_next(posfile) == 0);
    icheckEqual("last  block", 65600, posfile->current_block);

    printf("posfile tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
