/*  test/t_chrsplit.c -- chrsplit test cases.

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
#include <htslib/sam.h>

#include "bambi.h"

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

int main_chrsplit(int argc, char *argv[]);

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

int va_cmp(const void *a, const void *b)
{
    const char *s1 = *(char**)a;
    const char *s2 = *(char**)b;
    return strcmp(s1,s2);
}

static void va_sort(va_t *va)
{
    qsort(va->entries, va->end, sizeof(char*), va_cmp);
}

static bool va_compare(va_t *va1, va_t *va2)
{
    int n;
    if (va1->end != va2->end) return false;
    for (n=0; n<va1->end; n++) {
        if (strcmp((char*)(va1->entries[n]),(char*)(va2->entries[n])) != 0) return false;
    }
    return true;
}

static void va_dump(FILE *f, va_t *va)
{
    int n;
    for (n=0; n < va->end; n++) {
        fprintf(f,"\t%s\n", (char *)(va->entries[n]));
    }
}


static void checkReadNames(va_t *expected, char *fname)
{
    samFile *f = hts_open_format(fname, "r", NULL);
    bam_hdr_t *h = sam_hdr_read(f);
    bam1_t *rec = bam_init1();
    va_t *got = va_init(10,free);

    while (sam_read1(f,h,rec) != -1) {
        char *qname = bam_get_qname(rec);
        if (va_contains(got,qname)==-1) va_push(got,strdup(qname));
    }
    va_sort(expected); va_sort(got);
    if (!va_compare(expected,got)) {
        failure++;
        fprintf(stderr,"checkReadNames for file %s failed:\nExpected\n", fname);
        va_dump(stderr,expected);
        fprintf(stderr,"Got:\n");
        va_dump(stderr,got);
    }
    va_free(got);
    hts_close(f);
    bam_hdr_destroy(h);
    bam_destroy1(rec);
}

static void testxahuman(char *TMPDIR)
{
    int argc;
    char** argv;
    char target[512];
    char exclude[512];
    va_t *expected_target = va_init(10,NULL);
    va_t *expected_exclude = va_init(10,NULL);

    va_push(expected_target,"MT_MT");
    va_push(expected_target,"y_and_y");
    va_push(expected_target,"pair_unmapped");
    va_push(expected_target,"first_unmapped");
    va_push(expected_target,"second_unmapped");

    va_push(expected_exclude,"first_chimeric");
    va_push(expected_exclude,"twenty_twenty");
    va_push(expected_exclude,"unmapped_other");
    va_push(expected_exclude,"other_unmapped");
    va_push(expected_exclude,"second_chimeric");

    sprintf(target,"%s/chrsplit_target_1.bam", TMPDIR);
    sprintf(exclude,"%s/chrsplit_exclude_1.bam", TMPDIR);

    argc = 0;
    argv = (char**)calloc(sizeof(char*), 100);
    argv[argc++] = strdup("bambi");
    argv[argc++] = strdup("chrsplit");
    argv[argc++] = strdup("-i");
    argv[argc++] = strdup(MKNAME(DATA_DIR,"/10503_1_fix_mate.sam"));
    argv[argc++] = strdup("--input-fmt");
    argv[argc++] = strdup("sam");
    argv[argc++] = strdup("-o");
    argv[argc++] = strdup(target);
    argv[argc++] = strdup("-e");
    argv[argc++] = strdup(exclude);
    main_chrsplit(argc-1, argv+1);

    checkReadNames(expected_target, target);
    checkReadNames(expected_exclude, exclude);

    free_args(argv);
    va_free(expected_target);
    va_free(expected_exclude);
}

static void testxahuman_exclude_unaligned(char *TMPDIR)
{
    int argc;
    char** argv;
    char target[512];
    char exclude[512];
    va_t *expected_target = va_init(10,NULL);
    va_t *expected_exclude = va_init(10,NULL);

    va_push(expected_target,"MT_MT");
    va_push(expected_target,"y_and_y");

    va_push(expected_exclude,"pair_unmapped");
    va_push(expected_exclude,"first_unmapped");
    va_push(expected_exclude,"second_unmapped");
    va_push(expected_exclude,"first_chimeric");
    va_push(expected_exclude,"twenty_twenty");
    va_push(expected_exclude,"unmapped_other");
    va_push(expected_exclude,"other_unmapped");
    va_push(expected_exclude,"second_chimeric");

    sprintf(target,"%s/chrsplit_target_1.bam", TMPDIR);
    sprintf(exclude,"%s/chrsplit_exclude_1.bam", TMPDIR);

    argc = 0;
    argv = (char**)calloc(sizeof(char*), 100);
    argv[argc++] = strdup("bambi");
    argv[argc++] = strdup("chrsplit");
    argv[argc++] = strdup("-i");
    argv[argc++] = strdup(MKNAME(DATA_DIR,"/10503_1_fix_mate.sam"));
    argv[argc++] = strdup("--input-fmt");
    argv[argc++] = strdup("sam");
    argv[argc++] = strdup("-o");
    argv[argc++] = strdup(target);
    argv[argc++] = strdup("-e");
    argv[argc++] = strdup(exclude);
    argv[argc++] = strdup("-u");
    main_chrsplit(argc-1, argv+1);

    checkReadNames(expected_target, target);
    checkReadNames(expected_exclude, exclude);

    free_args(argv);
    va_free(expected_target);
    va_free(expected_exclude);
}

static void testyhuman(char *TMPDIR)
{
    int argc;
    char** argv;
    char target[512];
    char exclude[512];
    va_t *expected_target = va_init(10,NULL);
    va_t *expected_exclude = va_init(10,NULL);

    va_push(expected_target,"MT_MT");
    va_push(expected_target,"twenty_twenty");
    va_push(expected_target,"unmapped_other");
    va_push(expected_target,"pair_unmapped");
    va_push(expected_target,"other_unmapped");

    va_push(expected_exclude,"y_and_y");
    va_push(expected_exclude,"first_unmapped");
    va_push(expected_exclude,"second_unmapped");
    va_push(expected_exclude,"first_chimeric");
    va_push(expected_exclude,"second_chimeric");

    sprintf(target,"%s/chrsplit_target_1.bam", TMPDIR);
    sprintf(exclude,"%s/chrsplit_exclude_1.bam", TMPDIR);

    argc = 0;
    argv = (char**)calloc(sizeof(char*), 100);
    argv[argc++] = strdup("bambi");
    argv[argc++] = strdup("chrsplit");
    argv[argc++] = strdup("-i");
    argv[argc++] = strdup(MKNAME(DATA_DIR,"/10503_1.sam"));
    argv[argc++] = strdup("--input-fmt");
    argv[argc++] = strdup("sam");
    argv[argc++] = strdup("-o");
    argv[argc++] = strdup(target);
    argv[argc++] = strdup("-e");
    argv[argc++] = strdup(exclude);
    argv[argc++] = strdup("-V");
    argv[argc++] = strdup("--subset");
    argv[argc++] = strdup("Y");
    main_chrsplit(argc-1, argv+1);

    checkReadNames(expected_target, target);
    checkReadNames(expected_exclude, exclude);

    free_args(argv);
    va_free(expected_target);
    va_free(expected_exclude);
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

    testxahuman(TMPDIR);
    testxahuman_exclude_unaligned(TMPDIR);
    testyhuman(TMPDIR);

    printf("chrsplit tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
