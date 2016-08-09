/*  test/i2b/bclfile.c -- bclfile test cases.

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

#include "../../src/bclfile.c"
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
        fprintf(stderr, "%s: Expected: %s \tGot: %s\n", name, expected, actual);
        failure++;
    }
}

void checkEqual(char *name, char *expected, char *actual)
{
    if (actual == NULL) actual = "<null>";
    if (strcmp(expected, actual)) {
        fprintf(stderr, "%s: Expected: %s \tGot: %s\n", name, expected, actual);
        failure++;
    }
}

void ccheckEqual(char *name, char expected, char actual)
{
    if (expected != actual) {
        fprintf(stderr, "%s: Expected: '%c' \tGot: '%c'\n", name, expected, actual);
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
    int n;

    bclfile_t *bclfile;

    // BCL tests

    bclfile = bclfile_open(MKNAME(DATA_DIR,"/../i2b/110323_HS13_06000_B_B039WABXX/Data/Intensities/BaseCalls/L001/C1.1/s_1_1101.bcl"));
    if (bclfile->errmsg) {
        fprintf(stderr,"Error opening file: %s\n", bclfile->errmsg);
        failure++;
    }
    icheckEqual("Total clusters", 2609912, bclfile->total_clusters);
    icheckEqual("Current cluster", 0, bclfile->current_cluster);

    bclfile_next(bclfile);
    ccheckEqual("Base", 'N', bclfile->base);
    icheckEqual("Quality", 0, bclfile->quality);
    icheckEqual("current cluster", 1, bclfile->current_cluster);

    for (n=0; n<306; n++) {
        bclfile_next(bclfile);
    }
    ccheckEqual("307 Base", 'A', bclfile->base);
    icheckEqual("307 Quality", 30, bclfile->quality);
    icheckEqual("307 current cluster", 307, bclfile->current_cluster);
    icheckEqual("307 Total clusters", 2609912, bclfile->total_clusters);

    while (bclfile_next(bclfile) == 0);
    ccheckEqual("last Base", 'G', bclfile->base);
    icheckEqual("last Quality", 20, bclfile->quality);
    icheckEqual("last current cluster", 2609912, bclfile->current_cluster);
    icheckEqual("last Total clusters", 2609912, bclfile->total_clusters);

    bclfile_close(bclfile);

    // SCL tests

    bclfile = bclfile_open(MKNAME(DATA_DIR,"/../i2b/110323_HS13_06000_B_B039WABXX/Data/Intensities/BaseCalls/L001/C1.1/s_1_1101.scl"));
    if (bclfile->errmsg) {
        fprintf(stderr,"Error opening file: %s\n", bclfile->errmsg);
        failure++;
    }
    icheckEqual("SCL File Type", SCL, bclfile->file_type);
    icheckEqual("SCL Total clusters", 2609912, bclfile->total_clusters);
    icheckEqual("SCL Current cluster", 0, bclfile->current_cluster);

    bclfile_next(bclfile);
    ccheckEqual("SCL First Base", 'A', bclfile->base);

    for (n=0; n<306; n++) {
        bclfile_next(bclfile);
    }
    ccheckEqual("SCL 307 Base", 'T', bclfile->base);

    while (bclfile_next(bclfile) == 0);
    ccheckEqual("SCL Last Base", 'C', bclfile->base);

    printf("bclfile tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
