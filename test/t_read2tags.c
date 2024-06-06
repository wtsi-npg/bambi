/*  test/t_read2tags.c -- select test cases.

    Copyright (C) 2024 Genome Research Ltd.

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
#include "bamit.h"
#include <htslib/kstring.h>

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

int main_read2tags(int srgc, char *argv[]);

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
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("Ba");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("Qa");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("1:1:1");
}

void setup_test_2(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("Ba,Bb");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("Qa,Qb");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("1:2:4,1:3:5");
}

void setup_test_3(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("Ba");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("Qa");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("1:1:999");
    (*argv)[(*argc)++] = strdup("-d");
    (*argv)[(*argc)++] = strdup("ci");
}

void setup_test_4(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("Ba");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("Qa");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("2:1:999");
    (*argv)[(*argc)++] = strdup("-k");
    (*argv)[(*argc)++] = strdup("ci,RG");
}

void setup_test_5(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags_5.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("Ba");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("Qa");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("1:10");
}

void setup_test_6(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("Ba,Ba");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("Qa,Qb");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("1:2:2,1:1:1");
}

void setup_test_7(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("BC");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("QT");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("1:1:1");
    (*argv)[(*argc)++] = strdup("--replace");
}

void setup_test_8(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("BC");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("QT");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("1:1:1");
    (*argv)[(*argc)++] = strdup("--merge");
}

void setup_test_9(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("BC");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("QT");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("2:1:999");
    (*argv)[(*argc)++] = strdup("-d");
    (*argv)[(*argc)++] = strdup("ci,RG");
    (*argv)[(*argc)++] = strdup("-k");
    (*argv)[(*argc)++] = strdup("BC,QT");
    (*argv)[(*argc)++] = strdup("--merge");
}

void setup_test_10(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("BC");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("QT");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("2:1:999");
    (*argv)[(*argc)++] = strdup("-d");
    (*argv)[(*argc)++] = strdup("ci,RG");
    (*argv)[(*argc)++] = strdup("-k");
    (*argv)[(*argc)++] = strdup("BC,QT");
    (*argv)[(*argc)++] = strdup("--replace");
}

void setup_test_11(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("Ba");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("Qa");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("1:2:1:1");
}

void setup_test_12(int* argc, char*** argv, char *outputfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("read2tags");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-t");
    (*argv)[(*argc)++] = strdup("Ba");
    (*argv)[(*argc)++] = strdup("-q");
    (*argv)[(*argc)++] = strdup("Qa");
    (*argv)[(*argc)++] = strdup("-p");
    (*argv)[(*argc)++] = strdup("2:1:1:1");
}

void checkFiles(char *gotfile, char *expectfile, int verbose)
{
    BAMit_t *bgot = BAMit_open(gotfile, 'r', NULL, 0, NULL);
    BAMit_t *bexp = BAMit_open(expectfile, 'r', NULL, 0, NULL);
   // bam1_t *got_rec, *exp_rec;

    int f = failure;

    int c1 = sam_hdr_count_lines(bgot->h, "RG");
    int c2 = sam_hdr_count_lines(bexp->h, "RG");
    if (c1 != c2) {
        failure++; 
        if (verbose) fprintf(stderr, "RG lines: expected %d, got %d\n", c2, c1);
    }

    for (int n=0; n < c1; n++) {
        kstring_t ks_got; ks_initialize(&ks_got);
        kstring_t ks_exp; ks_initialize(&ks_exp);
        sam_hdr_find_line_pos(bgot->h, "RG", n, &ks_got);
        sam_hdr_find_line_pos(bexp->h, "RG", n, &ks_exp);
        if (strcmp(ks_str(&ks_got), ks_str(&ks_exp))) { 
            if (verbose) fprintf(stderr, "RG %d: expected %s, got %s\n", n, ks_str(&ks_exp), ks_str(&ks_got));
            failure++; 
            break; 
        }
        ks_free(&ks_got); ks_free(&ks_exp);
    }

    BAMit_free(bexp);
    BAMit_free(bgot);

    FILE *getfp = fopen(gotfile, "r");
    FILE *expfp = fopen(expectfile, "r");
    char getline[2048];
    char expline[2048];

    if (!getfp) {
        fprintf(stderr, "Can't open file %s\n", gotfile);
        exit(1);
    }

    if (!expfp) {
        fprintf(stderr, "Can't open file %s\n", expectfile);
        exit(1);
    }

    // skip header
    while (fgets(getline, 2047, getfp) > 0) {
        if (getline[0] != '@') break;
    }
    while (fgets(expline, 2047, expfp) > 0) {
        if (expline[0] != '@') break;
    }

    // compare read records
    while (true) {
        if (strcmp(getline,expline) != 0) {
            fprintf(stderr, "Expected: %sFound   : %s\n", expline, getline);
            failure++;
        }
        if (fgets(getline, 2047, getfp) == 0) break;
        if (fgets(expline, 2047, expfp) == 0) break;
    }

    fclose(getfp); fclose(expfp);

    if (verbose) {
        if (f == failure) fprintf(stderr, " :\tpass\n");
        else              fprintf(stderr, " :\t*** FAIL ***\n");
    }

    return;
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

    int argc_1;
    char** argv_1;
    char outputfile[1024];

    // minimal options
    if (verbose) fprintf(stderr,"Test 1: minimal options\n");
    sprintf(outputfile,"%s/read2tags_1.sam", TMPDIR);
    setup_test_1(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_1.sam"),verbose);
    free_args(argv_1);

    // overlapping reads
    if (verbose) fprintf(stderr,"Test 2: Overlapping reads\n");
    sprintf(outputfile,"%s/read2tags_2.sam", TMPDIR);
    setup_test_2(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_2.sam"),verbose);
    free_args(argv_1);

    // remove first record
    if (verbose) fprintf(stderr,"Test 3: remove first record\n");
    sprintf(outputfile,"%s/read2tags_3.sam", TMPDIR);
    setup_test_3(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_3.sam"),verbose);
    free_args(argv_1);

    // remove second record
    if (verbose) fprintf(stderr,"Test 4: remove second record\n");
    sprintf(outputfile,"%s/read2tags_4.sam", TMPDIR);
    setup_test_4(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_4.sam"),verbose);
    free_args(argv_1);

    // handle single reads
    if (verbose) fprintf(stderr,"Test 5: handle single reads\n");
    sprintf(outputfile,"%s/read2tags_5.sam", TMPDIR);
    setup_test_5(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_5.sam"),verbose);
    free_args(argv_1);

    // specify duplicate tags
    if (verbose) fprintf(stderr,"Test 6: specify duplicate tags\n");
    sprintf(outputfile,"%s/read2tags_6.sam", TMPDIR);
    setup_test_6(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_6.sam"),verbose);
    free_args(argv_1);

    // use --replace option
    if (verbose) fprintf(stderr,"Test 7: use --replace option\n");
    sprintf(outputfile,"%s/read2tags_7.sam", TMPDIR);
    setup_test_7(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_7.sam"),verbose);
    free_args(argv_1);

    // use --merge option
    if (verbose) fprintf(stderr,"Test 8: use --merge option\n");
    sprintf(outputfile,"%s/read2tags_8.sam", TMPDIR);
    setup_test_8(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_8.sam"),verbose);
    free_args(argv_1);

    // use --merge option with duplicate tags
    if (verbose) fprintf(stderr,"Test 9: use --merge option with duplicate tags\n");
    sprintf(outputfile,"%s/read2tags_9.sam", TMPDIR);
    setup_test_9(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_9.sam"),verbose);
    free_args(argv_1);

    // use --replace option with duplicate tags
    if (verbose) fprintf(stderr,"Test 10: use --replace option with duplicate tags\n");
    sprintf(outputfile,"%s/read2tags_10.sam", TMPDIR);
    setup_test_10(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_10.sam"),verbose);
    free_args(argv_1);

    // write tags to read 2 from read 1
    if (verbose) fprintf(stderr,"Test 11: write tags to read 2 from read 1\n");
    sprintf(outputfile,"%s/read2tags_11.sam", TMPDIR);
    setup_test_11(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_11.sam"),verbose);
    free_args(argv_1);

    // write tags to read 1 from read 2
    if (verbose) fprintf(stderr,"Test 12: write tags to read 1 from read 2\n");
    sprintf(outputfile,"%s/read2tags_12.sam", TMPDIR);
    setup_test_12(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_12.sam"),verbose);
    free_args(argv_1);

    printf("read2tags tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
