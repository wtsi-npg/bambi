/*  test/t_adapters.c -- select test cases.

    Copyright (C) 2019 Genome Research Ltd.

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
#include <unistd.h>
#include <stdio.h>
#include <math.h>

#include "bamit.h"
#include "bambi.h"

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

int main_adapters(int argc, char *argv[]);

const char * bambi_version(void)
{
    return "12.34";
}

int success = 0;
int failure = 0;

static void dumpTags1(bam1_t *rec, char *title)
{
    uint8_t *s;
    fprintf(stderr, "%s: %s\t", title, bam_get_qname(rec));
    s = bam_aux_get(rec,"aa");
    fprintf(stderr, "aa: %s\t", (s ? bam_aux2Z(s) : "NULL"));
    s = bam_aux_get(rec,"af");
    fprintf(stderr, "af: %f\t", (s ? bam_aux2f(s) : 0));
    s = bam_aux_get(rec,"ar");
    fprintf(stderr, "ar: %f\t", (s ? bam_aux2f(s) : 0));
    s = bam_aux_get(rec,"as");
    fprintf(stderr, "as: %d\t", (s ? (int)bam_aux2i(s) : 0));
    s = bam_aux_get(rec,"ah");
    fprintf(stderr, "ah: %d\t", (s ? (int)bam_aux2i(s) : 0));
    s = bam_aux_get(rec,"a3");
    fprintf(stderr, "a3: %d\t", (s ? (int)bam_aux2i(s) : 0));
    s = bam_aux_get(rec,"sc");
    if (s) fprintf(stderr, "sc: %d\t", (int)bam_aux2i(s));

    fprintf(stderr,"\n");
}

static void dumpTags(bam1_t *gotrec, bam1_t *exprec)
{
    dumpTags1(exprec, "Expected");
    dumpTags1(gotrec, "Found   ");
    fprintf(stderr,"\n");
}

void compareTags(bam1_t *gotrec, bam1_t *exprec, char *tag, char tagtype)
{
    uint8_t *s;
    char *gottag_Z=NULL, *exptag_Z=NULL;
    int gottag_i=0, exptag_i=0;
    double gottag_f=0, exptag_f=0;
    bool diff = false;

    s = bam_aux_get(gotrec,"aa"); if (s && *s=='Z') gottag_Z = bam_aux2Z(s);
    s = bam_aux_get(exprec,"aa"); if (s && *s=='Z') exptag_Z = bam_aux2Z(s);
    if ( !(gottag_Z==NULL && exptag_Z==NULL) ) {
        if ( (gottag_Z && !exptag_Z) ||
             (!gottag_Z && exptag_Z) ||
             (gottag_Z && exptag_Z && strcmp(gottag_Z,exptag_Z) != 0) ) { 
            diff = true; failure++;
        }
    }
    s = bam_aux_get(gotrec,"as"); if (s) gottag_i = bam_aux2i(s);
    s = bam_aux_get(exprec,"as"); if (s) exptag_i = bam_aux2i(s);
    if ( !(gottag_i==0 && exptag_i==0) ) {
        if ( (gottag_i && !exptag_i) ||
             (!gottag_i && exptag_i) ||
             (gottag_i && exptag_i && gottag_i != exptag_i) ) { 
            diff = true; failure++;
        }
    }

    s = bam_aux_get(gotrec,"ah"); if (s) gottag_i = bam_aux2i(s);
    s = bam_aux_get(exprec,"ah"); if (s) exptag_i = bam_aux2i(s);
    if ( !(gottag_i==0 && exptag_i==0) ) {
        if ( (gottag_i && !exptag_i) ||
             (!gottag_i && exptag_i) ||
             (gottag_i && exptag_i && gottag_i != exptag_i) ) { 
            diff = true; failure++;
        }
    }

    s = bam_aux_get(gotrec,"a3"); if (s) gottag_i = bam_aux2i(s);
    s = bam_aux_get(exprec,"a3"); if (s) exptag_i = bam_aux2i(s);
    if ( !(gottag_i==0 && exptag_i==0) ) {
        if ( (gottag_i && !exptag_i) ||
             (!gottag_i && exptag_i) ||
             (gottag_i && exptag_i && gottag_i != exptag_i) ) { 
            diff = true; failure++;
        }
    }

    s = bam_aux_get(gotrec,"af"); if (s && *s=='f') gottag_f = bam_aux2f(s);
    s = bam_aux_get(exprec,"af"); if (s && *s=='f') exptag_f = bam_aux2f(s);
    if ( !(gottag_f==0 && exptag_f==0) ) {
        if ( (gottag_f && !exptag_f) ||
             (!gottag_f && exptag_f) ||
             (gottag_f && exptag_f && fabs(gottag_f-exptag_f)>.0001 ) ) { 
            diff = true; failure++;
        }
    }

    s = bam_aux_get(gotrec,"ar"); if (s && *s=='f') gottag_f = bam_aux2f(s);
    s = bam_aux_get(exprec,"ar"); if (s && *s=='f') exptag_f = bam_aux2f(s);
    if ( !(gottag_f==0 && exptag_f==0) ) {
        if ( (gottag_f && !exptag_f) ||
             (!gottag_f && exptag_f) ||
             (gottag_f && exptag_f && fabs(gottag_f-exptag_f)>.0001 ) ) { 
            diff = true; failure++;
        }
    }
    if (diff) dumpTags(gotrec, exprec);
}

void checkFiles(char *gotfile, char *expectfile, int verbose)
{
    BAMit_t *bgot = BAMit_open(gotfile, 'r', NULL, 0, NULL);
    BAMit_t *bexp = BAMit_open(expectfile, 'r', NULL, 0, NULL);

    bam1_t *got_rec, *exp_rec;

    while ((exp_rec = BAMit_next(bexp)) != NULL) {
        got_rec = BAMit_next(bgot);
        if (!got_rec) { fprintf(stderr, "%s ended too soon\n", gotfile); failure++; return; }
        // check qname
        if (strcmp(bam_get_qname(exp_rec), bam_get_qname(exp_rec)) != 0) {
            fprintf(stderr,"Qname differs: expected: %s\n", bam_get_qname(exp_rec));
            fprintf(stderr,"               got     : %s\n", bam_get_qname(got_rec));
            failure++;
            break;
        }

        // check tags
        compareTags(got_rec, exp_rec, "aa", 'Z');
    }

    BAMit_free(bexp);
    BAMit_free(bgot);
    return;
}

void test_1(char *TMPDIR, int verbose)
{
    char outputfile[1024];
    char cmd[1024];

    sprintf(outputfile,"%s/adapters_1.bam", TMPDIR);
    sprintf(cmd, "src/bambi adapters -o %s/adapters_1.bam %s", TMPDIR, MKNAME(DATA_DIR,"/adapters.bam"));
    if (system(cmd)) { fprintf(stderr, "Command %s failed\n", cmd); failure++; return; }
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/adapters.bam"),verbose);
}

int main(int argc, char**argv)
{
    int verbose = 0;

    int getopt_char;
    while ((getopt_char = getopt(argc, argv, "v")) != -1) {
        switch (getopt_char) {
            case 'v': ++verbose;
                      break;
            default: printf("usage: t_adapters [-v]\n\n"
                            " -v verbose output\n"
                           );
                     break;
        }
    }

    // Cleanup getopt
    optind = 1;

    // create temp directory
    char template[] = "/tmp/bambi.XXXXXX";
    char *TMPDIR = mkdtemp(template);
    if (TMPDIR == NULL) {
        fprintf(stderr,"Can't create temp directory\n");
        exit(1);
    } else {
        if (verbose) fprintf(stderr,"Created temporary directory: %s\n", TMPDIR);
    }

    test_1(TMPDIR, verbose);

    printf("adapters tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
