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

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "bclfile.h"

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

void display(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}

void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap,fmt);
    fflush(stdout);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

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

    bclfile = bclfile_open(MKNAME(DATA_DIR,"/s_1_1101.bcl.gz"), MT_HISEQX, -1);
    if (bclfile->errmsg) {
        fprintf(stderr,"Error opening file: %s\n", bclfile->errmsg);
        failure++;
    }
    icheckEqual("Total clusters", 2609912, bclfile->total_clusters);
    icheckEqual("Current cluster", 0, bclfile->base_ptr);

    ccheckEqual("Base", 'N', bclfile_base(bclfile,0));
    icheckEqual("Quality", 0, bclfile_quality(bclfile,0));

    ccheckEqual("307 Base", 'A', bclfile_base(bclfile,306));
    icheckEqual("307 Quality", 30, bclfile_quality(bclfile,306));

    bclfile_close(bclfile);

    // CBCL tests

    bclfile = bclfile_open(MKNAME(DATA_DIR,"/novaseq/Data/Intensities/BaseCalls/L001/C1.1/L001_1.cbcl"), MT_NOVASEQ, 1101);
    if (bclfile->errmsg) {
        fprintf(stderr,"Error opening file: %s\n", bclfile->errmsg);
        failure++;
    }
    bclfile_load_tile(bclfile,1101,NULL,-1);
    icheckEqual("machine type", MT_NOVASEQ, bclfile->machine_type);
    icheckEqual("version", 1, bclfile->version);
    icheckEqual("header_size", 65, bclfile->header_size);
    icheckEqual("bits_per_base", 2, bclfile->bits_per_base);
    icheckEqual("bits_per_qual", 2, bclfile->bits_per_qual);
    icheckEqual("number of bins", 4, bclfile->nbins);
    icheckEqual("number of tiles", 1, bclfile->ntiles);

    for (n=0; n < bclfile->ntiles; n++) {
        tilerec_t *t = bclfile->tiles->entries[n];
        fprintf(stderr,"%d\t%d\t%d\t%d\t%d\n", t->tilenum, t->nclusters, t->uncompressed_blocksize, t->compressed_blocksize, bclfile->pfFlag);
    }

    ccheckEqual("CBCL First Base", 'T', bclfile_base(bclfile,0));
    ccheckEqual("CBCL Second Base", 'G', bclfile_base(bclfile,1));
    ccheckEqual("CBCL Third Base", 'N', bclfile_base(bclfile,2));
    ccheckEqual("CBCL Last Base", 'G', bclfile_base(bclfile,27));

    icheckEqual("CBCL Number of bases", 28, bclfile->bases_size);

    bclfile_close(bclfile);

    printf("bclfile tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
