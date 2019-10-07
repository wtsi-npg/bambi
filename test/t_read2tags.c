/*  test/t_select.c -- select test cases.

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
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags_5.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/read2tags.sam"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    sprintf(outputfile,"%s/read2tags_1.bam", TMPDIR);
    setup_test_1(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_1.bam"),verbose);
    free_args(argv_1);

    // overlapping reads
    sprintf(outputfile,"%s/read2tags_2.bam", TMPDIR);
    setup_test_2(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_2.bam"),verbose);
    free_args(argv_1);

    // remove first record
    sprintf(outputfile,"%s/read2tags_3.bam", TMPDIR);
    setup_test_3(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_3.bam"),verbose);
    free_args(argv_1);

    // remove second record
    sprintf(outputfile,"%s/read2tags_4.bam", TMPDIR);
    setup_test_4(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_4.bam"),verbose);
    free_args(argv_1);

    // handle single reads
    sprintf(outputfile,"%s/read2tags_5.bam", TMPDIR);
    setup_test_5(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_5.bam"),verbose);
    free_args(argv_1);

    // specify duplicate tags
    sprintf(outputfile,"%s/read2tags_6.bam", TMPDIR);
    setup_test_6(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_6.bam"),verbose);
    free_args(argv_1);

    // use --replace option
    sprintf(outputfile,"%s/read2tags_7.bam", TMPDIR);
    setup_test_7(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_7.bam"),verbose);
    free_args(argv_1);

    // use --merge option
    sprintf(outputfile,"%s/read2tags_8.bam", TMPDIR);
    setup_test_8(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_8.bam"),verbose);
    free_args(argv_1);

    // use --merge option with duplicate tags
    sprintf(outputfile,"%s/read2tags_9.bam", TMPDIR);
    setup_test_9(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_9.bam"),verbose);
    free_args(argv_1);

    // use --replace option with duplicate tags
    sprintf(outputfile,"%s/read2tags_10.bam", TMPDIR);
    setup_test_10(&argc_1, &argv_1, outputfile);
    main_read2tags(argc_1-1, argv_1+1);
    checkFiles(outputfile,MKNAME(DATA_DIR,"/out/read2tags_10.bam"),verbose);
    free_args(argv_1);


    printf("read2tags tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
