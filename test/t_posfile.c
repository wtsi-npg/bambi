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

#include "../src/posfile.c"
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
    int n;

    posfile_t *posfile;

    /*
     * clocs file
     */
//    posfile = posfile_open(MKNAME(DATA_DIR,"/s_1_1101.clocs"));
    posfile = posfile_open(MKNAME(DATA_DIR,"/test.clocs"));
    if (posfile->errmsg) {
        fprintf(stderr,"Error opening file: %s\n", posfile->errmsg);
        failure++;
    }
    posfile_load(posfile,0,NULL);

    icheckEqual("Version", 1, posfile->version);
    icheckEqual("Total blocks", 65600, posfile->total_blocks);

    icheckEqual("next X", 1235, posfile_get_x(posfile,0));
    icheckEqual("next Y", 1989, posfile_get_y(posfile,0));

    icheckEqual("307 x", 1279, posfile_get_x(posfile,306));
    icheckEqual("307 y", 2120, posfile_get_y(posfile,306));

    posfile_close(posfile);

    /*
     * locs file
     */
    posfile = posfile_open(MKNAME(DATA_DIR,"/test.locs"));
    if (posfile->errmsg) {
        fprintf(stderr,"Error opening file: %s\n", posfile->errmsg);
        failure++;
    }
    posfile_load(posfile,0,NULL);

    icheckEqual("LOCS: Total blocks", 500, posfile->total_blocks);
    icheckEqual("LOCS: current block", 0, posfile->current_block);

    icheckEqual("LOCS: first X", 6568, posfile_get_x(posfile,0));
    icheckEqual("LOCS: first Y", 21421, posfile_get_y(posfile,0));

    icheckEqual("LOCS: last x", 19845, posfile_get_x(posfile,499));
    icheckEqual("LOCS: last y", 7503, posfile_get_y(posfile,499));

    posfile_close(posfile);

    printf("posfile tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
