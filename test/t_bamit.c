/*  t_bamit.c -- bamit test cases.

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
#include "bamit.h"
#include "htslib/hts.h"

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

void dumpRec(bam1_t *rec)
{
    fprintf(stderr,"mpos: %d  pos: %d  qname: %s  flag: %d\n", rec->core.mpos, rec->core.pos, bam_get_qname(rec), rec->core.flag);
}

int main(int argc, char**argv)
{
    int n;

    hts_verbose=6;

    bam1_t *rec;
    samFile *f = hts_open_format(MKNAME(DATA_DIR,"/bamit.bam"), "rb", NULL);
    bam_hdr_t *h = sam_hdr_read(f);

    BAMit_t *bit = BAMit_init(f,h);
    icheckEqual("First record has next", true, BAMit_hasnext(bit));
    rec = BAMit_next(bit);
    checkEqual("First name", "IL16_986:1:9:9:307", bam_get_qname(rec));
    icheckEqual("First flag", 83, rec->core.flag);
    rec = BAMit_next(bit);
    checkEqual("Second name", "IL16_986:1:9:9:307", bam_get_qname(rec));
    icheckEqual("Second flag", 163, rec->core.flag);
    rec = BAMit_peek(bit);
    checkEqual("Peek name", "IL16_986:1:9:9:47", bam_get_qname(rec));
    icheckEqual("Peek flag", 83, rec->core.flag);
    rec = BAMit_next(bit);
    checkEqual("Third name", "IL16_986:1:9:9:47", bam_get_qname(rec));
    icheckEqual("Third flag", 83, rec->core.flag);
    
    n=2;
    while (rec) {
        n++; rec = BAMit_next(bit);
    }
    icheckEqual("Number of records", 6, n);
    icheckEqual("End of records", false, BAMit_hasnext(bit));
    BAMit_free(bit);

    bit = BAMit_open(MKNAME(DATA_DIR,"/bamit.bam"), 'r', NULL, 0);
    rec = BAMit_next(bit);
    checkEqual("First name", "IL16_986:1:9:9:307", bam_get_qname(rec));
    BAMit_free(bit);

    printf("BAMit tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
