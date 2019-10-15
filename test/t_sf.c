/*  test/t_sf.c -- spatial_filter test cases.

    Copyright (C) 2017 Genome Research Ltd.

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
#include <string.h>

#include "bamit.h"
#include <htslib/kstring.h>

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

int main_spatial_filter(int argc, char *argv[]);

const char * bambi_version(void)
{
    return "12.34";
}

int success = 0;
int failure = 0;

void checkFiles(char *gotfile, char *expectfile, int verbose)
{
    BAMit_t *bgot = BAMit_open(gotfile, 'r', NULL, 0, NULL);
    BAMit_t *bexp = BAMit_open(expectfile, 'r', NULL, 0, NULL);
    bam1_t *got_rec, *exp_rec;

    int c = sam_hdr_count_lines(bgot->h, "RG");
    if (c != sam_hdr_count_lines(bexp->h, "RG")) { failure++; return; }

    for (int n=0; n < c; n++) {
        kstring_t ks_got; ks_initialize(&ks_got);
        kstring_t ks_exp; ks_initialize(&ks_exp);
        sam_hdr_find_line_pos(bgot->h, "RG", n, &ks_got);
        sam_hdr_find_line_pos(bexp->h, "RG", n, &ks_exp);
        if (strcmp(ks_str(&ks_got), ks_str(&ks_exp))) { failure++; return; }
        ks_free(&ks_got); ks_free(&ks_exp);
    }

    while ((exp_rec = BAMit_next(bexp)) != NULL) {
        got_rec = BAMit_next(bgot);
        if (!got_rec) { fprintf(stderr, "%s ended too soon\n", gotfile); failure++; return; }
        if (memcmp(got_rec->data, exp_rec->data, got_rec->l_data)) {
            failure++;
            break;
        }
    }

    BAMit_free(bexp);
    BAMit_free(bgot);
    return;
}

void checkFilterFiles(char *prog, char *tmpdir, char *gotfile, char *expectfile)
{
    char cmd[1024];
    sprintf(cmd,"%s -D -v -F %s > %s/got.txt", prog, gotfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"%s -D -v -F %s > %s/expect.txt", prog, expectfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"diff -ICommand %s/got.txt %s/expect.txt", tmpdir, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
}

int main(int argc, char**argv)
{
    int verbose = 0;

    int getopt_char;
    while ((getopt_char = getopt(argc, argv, "v")) != -1) {
        switch (getopt_char) {
            case 'v': ++verbose;
                      break;
            default: printf("usage: test_parse_args [-v]\n\n"
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

    // minimal options
    char outputfile[512];
    char filterfile[512];
    char cmd[2048];
    char prog[512];

    sprintf(prog, "%s", "src/bambi spatial_filter");

    // create filter
    if (verbose) fprintf(stderr,"Creating filter\n");
    snprintf(filterfile, sizeof(filterfile), "%s/sf_1.filter", TMPDIR);
    snprintf(outputfile, sizeof(outputfile), "%s/sf_filtered.bam", TMPDIR);
    snprintf(cmd, sizeof(cmd), "%s -c -F %s %s", prog, filterfile, MKNAME(DATA_DIR,"/sf.bam"));
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    else checkFilterFiles(prog, TMPDIR, filterfile, MKNAME(DATA_DIR,"/out/sf_1.filter"));

    // apply filter
    if (verbose) fprintf(stderr,"Applying filter\n");
    snprintf(cmd, sizeof(cmd), "%s -a --verbose -F %s -o %s %s", prog, filterfile, outputfile, MKNAME(DATA_DIR,"/sf.bam"));
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    else checkFiles(outputfile, MKNAME(DATA_DIR,"/out/sf_filtered.bam"), verbose);

    printf("spatial_filter tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
