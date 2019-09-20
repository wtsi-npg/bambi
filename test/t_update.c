/*  test/t_update.c -- update test cases.

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
#include <unistd.h>
#include <stdio.h>

#include "bambi.h"

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

int main_update(int srgc, char *argv[]);

const char * bambi_version(void)
{
    return "12.34";
}

int success = 0;
int failure = 0;

void free_args(char **argv)
{
    for (int n=0; n<100; n++) free(argv[n]);
    free(argv);
}

void setup_test_1(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("update");
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/update_1.sam"));
    (*argv)[(*argc)++] = strdup(outputfile);
}

void checkFiles(char *gotfile, char *expectfile, int verbose)
{
    char cmd[1024];

    if (verbose) fprintf(stderr,"\nComparing files: %s with %s\n", gotfile, expectfile);
    sprintf(cmd,"diff -I ID:bambi %s %s", gotfile, expectfile);
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
            default: printf("usage: t_update [-v]\n\n"
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

    int argc_1;
    char** argv_1;
    char outputfile[1024];

    // minimal options
    sprintf(outputfile,"%s/update_1.sam", TMPDIR);
    setup_test_1(&argc_1, &argv_1, outputfile);
    main_update(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/update_1.sam"),verbose);
    free_args(argv_1);

    printf("read2tags tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
