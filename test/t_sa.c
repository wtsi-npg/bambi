/*  test/t_sa.c -- substitution_anlysis test cases.

    Copyright (C) 2021 Genome Research Ltd.

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

int main_substitution_anlysis(int argc, char *argv[]);

const char * bambi_version(void)
{
    return "12.34";
}

int success = 0;
int failure = 0;

void checkOutputFiles(char *tmpdir, char *gotfile, char *expectfile)
{
    char cmd[1024];
    sprintf(cmd,"diff %s/%s %s", tmpdir, gotfile, expectfile);
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
    char cmd[2048];
    char prog[512];

    sprintf(prog, "%s", "src/bambi substitution_analysis");

    snprintf(cmd, sizeof(cmd), "%s -o %s/sa.txt %s", prog, TMPDIR, MKNAME(DATA_DIR,"/sa.bam"));
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    else checkOutputFiles(TMPDIR, "sa.txt", MKNAME(DATA_DIR,"/out/sa.txt"));

    printf("substitution_analysis tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
