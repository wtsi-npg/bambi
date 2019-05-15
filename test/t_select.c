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

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

int main_select(int argc, char *argv[]);

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


void setup_test_1(int* argc, char*** argv, char *outputfile, char *metricsfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/select_1.sam"));
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/select_1_human.sam"));
    (*argv)[(*argc)++] = strdup("--input-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("-m");
    (*argv)[(*argc)++] = strdup(metricsfile);
}

void setup_test_2(int* argc, char*** argv, char *outputfile, char *unalignedfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/select_1.sam"));
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/select_1_human_unmapped_with_ref.sam"));
    (*argv)[(*argc)++] = strdup("--input-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("-n");
    (*argv)[(*argc)++] = strdup(unalignedfile);
}

void setup_test_3(int* argc, char*** argv, char *outputfile, char *unalignedfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/select_single.sam"));
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/select_single_human_unmapped_with_ref.sam"));
    (*argv)[(*argc)++] = strdup("--input-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("-n");
    (*argv)[(*argc)++] = strdup(unalignedfile);
}

void setup_test_4(int* argc, char*** argv, char *outputfile, char *unalignedfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/select_single.sam"));
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/select_single_human_with_sup.sam"));
    (*argv)[(*argc)++] = strdup("--input-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("-n");
    (*argv)[(*argc)++] = strdup(unalignedfile);
}

void setup_test_5(int* argc, char*** argv, char *outputfile, char *unalignedfile, char *metricsfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("select");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/chimeric.sam"));
    (*argv)[(*argc)++] = strdup("--input-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("-n");
    (*argv)[(*argc)++] = strdup(unalignedfile);
    (*argv)[(*argc)++] = strdup("-m");
    (*argv)[(*argc)++] = strdup(metricsfile);
}

void checkFiles(char *tmpdir, char *gotfile, char *expectfile, int verbose)
{
    char cmd[1024];

    if (verbose) fprintf(stderr,"\nComparing headers: %s with %s\n", gotfile, expectfile);
    // compare headers
    sprintf(cmd,"samtools view -H %s/%s |grep -v ^@PG |sort | perl -n -e 'chomp; @x=split /\t/;@y=sort @x; print join \",\",@y; print \"\n\";' | sed s:/tmp/bambi[^/]*:/tmp/xyzzy:g > %s/got.txt", tmpdir, gotfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"samtools view -H %s | grep -v ^@PG| sort | perl -n -e 'chomp; @x=split /\t/;@y=sort @x; print join \",\",@y; print \"\n\";' | sed s:/tmp/bambi[^/]*:/tmp/xyzzy:g > %s/expect.txt", expectfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"diff %s/got.txt %s/expect.txt", tmpdir, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }

    // compare records
    if (verbose) fprintf(stderr,"\nComparing records: %s with %s\n", gotfile, expectfile);
    sprintf(cmd,"samtools view %s/%s > %s/got.txt", tmpdir, gotfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"samtools view %s > %s/expect.txt", expectfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"diff %s/got.txt %s/expect.txt", tmpdir, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
}

void checkJSONFiles(char *tmpdir, char *gotfile, char *expectfile)
{
    char cmd[1024];
    sprintf(cmd,"cat %s | sed s/[,\\\\[\\\\{\\\\}]/\\\\n/g | grep -v tmp/bambi | sort > %s/got.txt", gotfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command filed: %s\n", cmd); failure++; }
    sprintf(cmd,"cat %s | sed s/[,\\\\[\\\\{\\\\}]/\\\\n/g | grep -v tmp/bambi | sort > %s/expect.txt", expectfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command filed: %s\n", cmd); failure++; }
    sprintf(cmd,"diff -B %s/got.txt %s/expect.txt", tmpdir, tmpdir);
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
    int argc_1;
    char** argv_1;
    char outputfile[1024];
    char unalignedfile[512];
    char metricsfile[512];

    sprintf(outputfile,"%s/select_1.bam,%s/select_1_human.bam", TMPDIR, TMPDIR);
    sprintf(metricsfile,"%s/select_1_metrics.json", TMPDIR);

    setup_test_1(&argc_1, &argv_1, outputfile, metricsfile);
    main_select(argc_1-1, argv_1+1);

    checkFiles(TMPDIR,"select_1.bam",MKNAME(DATA_DIR,"/out/select_1.bam"),verbose);
    checkFiles(TMPDIR,"select_1_human.bam",MKNAME(DATA_DIR,"/out/select_1_human.bam"),verbose);

    checkJSONFiles(TMPDIR,metricsfile,MKNAME(DATA_DIR,"/out/select_1_metrics.json"));
    free_args(argv_1);

    // unaligned file test

    sprintf(outputfile,"%s/select_2.bam,%s/select_2_human.bam", TMPDIR, TMPDIR);
    sprintf(unalignedfile,"%s/select_2_unaligned.bam",TMPDIR);

    setup_test_2(&argc_1, &argv_1, outputfile, unalignedfile);
    main_select(argc_1-1, argv_1+1);

    checkFiles(TMPDIR,"select_2.bam",MKNAME(DATA_DIR,"/out/select_2.bam"),verbose);
    checkFiles(TMPDIR,"select_2_human.bam",MKNAME(DATA_DIR,"/out/select_2_human.bam"),verbose);
    checkFiles(TMPDIR,"select_2_unaligned.bam",MKNAME(DATA_DIR,"/out/select_2_unaligned.bam"),verbose);
    free_args(argv_1);

    // single read data test

    sprintf(outputfile,"%s/select_3.bam,%s/select_3_human.bam", TMPDIR, TMPDIR);
    sprintf(unalignedfile,"%s/select_3_unaligned.bam",TMPDIR);

    setup_test_3(&argc_1, &argv_1, outputfile, unalignedfile);
    main_select(argc_1-1, argv_1+1);

    checkFiles(TMPDIR,"select_3.bam",MKNAME(DATA_DIR,"/out/select_3.bam"),verbose);
    checkFiles(TMPDIR,"select_3_human.bam",MKNAME(DATA_DIR,"/out/select_3_human.bam"),verbose);
    checkFiles(TMPDIR,"select_3_unaligned.bam",MKNAME(DATA_DIR,"/out/select_3_unaligned.bam"),verbose);
    free_args(argv_1);


    // supplemental read data test

    sprintf(outputfile,"%s/select_4.bam,%s/select_4_human.bam", TMPDIR, TMPDIR);
    sprintf(unalignedfile,"%s/select_4_unaligned.bam",TMPDIR);

    setup_test_4(&argc_1, &argv_1, outputfile, unalignedfile);
    main_select(argc_1-1, argv_1+1);

    checkFiles(TMPDIR,"select_4.bam",MKNAME(DATA_DIR,"/out/select_sup.bam"),verbose);
    checkFiles(TMPDIR,"select_4_human.bam",MKNAME(DATA_DIR,"/out/select_sup_human.bam"),verbose);
    checkFiles(TMPDIR,"select_4_unaligned.bam",MKNAME(DATA_DIR,"/out/select_sup_unaligned.bam"),verbose);
    free_args(argv_1);


    // chimeric metrics test

    sprintf(outputfile,"%s/select_5.bam", TMPDIR);
    sprintf(unalignedfile,"%s/select_5_unaligned.bam",TMPDIR);
    sprintf(metricsfile,"%s/select_5.json",TMPDIR);

    setup_test_5(&argc_1, &argv_1, outputfile, unalignedfile, metricsfile);
    main_select(argc_1-1, argv_1+1);

    checkJSONFiles(TMPDIR,metricsfile, MKNAME(DATA_DIR,"/out/chimeric.json"));
    free_args(argv_1);

    printf("select tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
