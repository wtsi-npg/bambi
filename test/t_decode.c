/*  test/decode/decode.c -- decode test cases.

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
#include "../src/hts_addendum.c"
#include "../src/decode.c"

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

#define NTHREADS 4

const char * bambi_version(void)
{
    return "12.34";
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

int success = 0;
int failure = 0;

static char *itoa(int i)
{
    const size_t sz = 32;
    char *a = malloc(sz);
    if (!a) die("Out of memory");
    snprintf(a, sz, "%d", i);
    return a;
}

void setup_test_1(int* argc, char*** argv, char *outputfile, char *metricsfile,
                  int threads)
{
    *argc = 16 + (threads ? 2 : 0);
    *argv = (char**)calloc(sizeof(char*), *argc);
    (*argv)[0] = strdup("bambi");
    (*argv)[1] = strdup("decode");
    (*argv)[2] = strdup("-i");
    (*argv)[3] = strdup(MKNAME(DATA_DIR,"/decode_1.sam"));
    (*argv)[4] = strdup("-o");
    (*argv)[5] = strdup(outputfile);
    (*argv)[6] = strdup("--output-fmt");
    (*argv)[7] = strdup("sam");
    (*argv)[8] = strdup("--input-fmt");
    (*argv)[9] = strdup("sam");
    (*argv)[10] = strdup("--barcode-file");
    (*argv)[11] = strdup(MKNAME(DATA_DIR,"/decode_1.tag"));
    (*argv)[12] = strdup("--metrics-file");
    (*argv)[13] = strdup(metricsfile);
    (*argv)[14] = strdup("--barcode-tag-name");
    (*argv)[15] = strdup("RT");
    if (threads) {
        (*argv)[16] = strdup("-t");
        (*argv)[17] = itoa(threads);
    }
}

void setup_test_2(int* argc, char*** argv, char *outputfile, char *metricsfile, char *chksumfile,
                  int threads)
{
    *argc = 20 + (threads ? 2 : 0);
    *argv = (char**)calloc(sizeof(char*), *argc);
    (*argv)[0] = strdup("bambi");
    (*argv)[1] = strdup("decode");
    (*argv)[2] = strdup("-i");
    (*argv)[3] = strdup(MKNAME(DATA_DIR,"/decode_1.sam"));
    (*argv)[4] = strdup("-o");
    (*argv)[5] = strdup(outputfile);
    (*argv)[6] = strdup("--output-fmt");
    (*argv)[7] = strdup("sam");
    (*argv)[8] = strdup("--input-fmt");
    (*argv)[9] = strdup("sam");
    (*argv)[10] = strdup("--barcode-file");
    (*argv)[11] = strdup(MKNAME(DATA_DIR,"/decode_1.tag"));
    (*argv)[12] = strdup("--convert-low-quality");
    (*argv)[13] = strdup("--change-read-name");
    (*argv)[14] = strdup("--metrics-file");
    (*argv)[15] = strdup(metricsfile);
    (*argv)[16] = strdup("--barcode-tag-name");
    (*argv)[17] = strdup("RT");
    (*argv)[18] = strdup("--chksum-file");
    (*argv)[19] = strdup(chksumfile);
    if (threads) {
        (*argv)[20] = strdup("-t");
        (*argv)[21] = itoa(threads);
    }
}

void setup_test_3(int* argc, char*** argv, char *outputfile, char *metricsfile, char *chksumfile,
                  int threads)
{
    *argc = 19 + (threads ? 2 : 0);
    *argv = (char**)calloc(sizeof(char*), *argc);
    (*argv)[0] = strdup("bambi");
    (*argv)[1] = strdup("decode");
    (*argv)[2] = strdup("-i");
    (*argv)[3] = strdup(MKNAME(DATA_DIR,"/decode_3.sam"));
    (*argv)[4] = strdup("-o");
    (*argv)[5] = strdup(outputfile);
    (*argv)[6] = strdup("--output-fmt");
    (*argv)[7] = strdup("sam");
    (*argv)[8] = strdup("--input-fmt");
    (*argv)[9] = strdup("sam");
    (*argv)[10] = strdup("--barcode-file");
    (*argv)[11] = strdup(MKNAME(DATA_DIR,"/decode_3.tag"));
    (*argv)[12] = strdup("--convert-low-quality");
    (*argv)[13] = strdup("--max-no-calls");
    (*argv)[14] = strdup("6");
    (*argv)[15] = strdup("--hash");
    (*argv)[16] = strdup("crc32");
    (*argv)[17] = strdup("--chksum-file");
    (*argv)[18] = strdup(chksumfile);
    if (threads) {
        (*argv)[19] = strdup("--threads");
        (*argv)[20] = itoa(threads);
    }
}

void setup_test_4(int* argc, char*** argv, char *outputfile, char* metricsfile,
                  int threads)
{
    *argc = 15 + (threads ? 2 : 0);
    *argv = (char**)calloc(sizeof(char*), *argc);
    (*argv)[0] = strdup("bambi");
    (*argv)[1] = strdup("decode");
    (*argv)[2] = strdup("-i");
    (*argv)[3] = strdup(MKNAME(DATA_DIR,"/decode_4.sam"));
    (*argv)[4] = strdup("-o");
    (*argv)[5] = strdup(outputfile);
    (*argv)[6] = strdup("--output-fmt");
    (*argv)[7] = strdup("sam");
    (*argv)[8] = strdup("--input-fmt");
    (*argv)[9] = strdup("sam");
    (*argv)[10] = strdup("--barcode-file");
    (*argv)[11] = strdup(MKNAME(DATA_DIR,"/decode_4.tag"));
    (*argv)[12] = strdup("--metrics-file");
    (*argv)[13] = strdup(metricsfile);
    (*argv)[14] = strdup("--ignore-pf");
    if (threads) {
        (*argv)[15] = strdup("--threads");
        (*argv)[16] = itoa(threads);
    }
}

void free_argv(int argc, char *argv[])
{
    for (int n=0; n < argc; free(argv[n++]));
    free(argv);
}

void test_noCalls(char *s, int e)
{
    int n;
    if ((n=noCalls(s)) == e) success++;
    else { failure++; fprintf(stderr, "noCalls(%s) returned %d: expected %d\n",s,n,e); }
}

void test_countMismatches(char *a, char *b, int e)
{
    int n;
    if ((n=countMismatches(a,b,999)) == e) success++;
    else { failure++; fprintf(stderr, "countMismatches(%s,%s) returned %d: expected %d\n", a,b,n,e); }
}

int main(int argc, char**argv)
{
    // test state
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

    // test isNoCall()
    if (!isNoCall('A')) success++;
    else { failure++; fprintf(stderr, "isNoCall('A') returned True\n"); }
    if (isNoCall('N')) success++;
    else { failure++; fprintf(stderr, "isNoCall('N') returned False\n"); }
    if (isNoCall('n')) success++;
    else { failure++; fprintf(stderr, "isNoCall('n') returned False\n"); }
    if (isNoCall('.')) success++;
    else { failure++; fprintf(stderr, "isNoCall('.') returned False\n"); }

    // test noCalls()
    test_noCalls("ABC",0);
    test_noCalls("ABCN",1);
    test_noCalls("NABCN",2);
    test_noCalls("NA+CN",2);

    // test countMismatches()
    test_countMismatches("ABC","AXC",1);
    test_countMismatches("ABC","XYZ",3);
    test_countMismatches("ABC","ABC",0);
    test_countMismatches("ABCxXYZ","ABCNXYZ",0);
    test_countMismatches("ABCiXYZ","ABCNXYZ",0);
    test_countMismatches("xBCiXYZ","NBCNXYz",1);
    test_countMismatches("AGCACGTT","AxCACGTTXXXXXX",1);

    //
    // Now test the actual decoding
    //

    unsigned int max_path_length = strlen(TMPDIR) + 100;
    char *outputfile = calloc(1,max_path_length);
    char *metricsfile = calloc(1,max_path_length);
    char *chksumfile = calloc(1,max_path_length);
    char cmd[1024];
    
    // minimal options

    for (int threads = 0; threads <= NTHREADS; threads+=NTHREADS) {
        int argc_1;
        char** argv_1;
        int result;

        snprintf(outputfile, max_path_length, "%s/decode_1%s.sam", TMPDIR, threads ? "threads" : "");
        snprintf(metricsfile, max_path_length, "%s/decode_1%s.metrics", TMPDIR, threads ? "threads" : "");
        setup_test_1(&argc_1, &argv_1, outputfile, metricsfile, threads);
        main_decode(argc_1-1, argv_1+1);
        free_argv(argc_1,argv_1);

        snprintf(cmd, sizeof(cmd), "diff -I ID:bambi %s %s", outputfile, MKNAME(DATA_DIR,"/out/6383_9_nosplit_nochange.sam"));
        result = system(cmd);
        if (result) {
            fprintf(stderr, "test 1 failed\n");
            failure++;
        } else {
            success++;
        }

        snprintf(cmd, sizeof(cmd), "diff -I ID:bambi %s %s", metricsfile, MKNAME(DATA_DIR,"/out/decode_1.metrics"));
        result = system(cmd);
        if (result) {
            fprintf(stderr, "test 1 failed at metrics file diff\n");
            failure++;
        } else {
            success++;
        }
    }

    // --convert_low_quality option
    for (int threads = 0; threads <= NTHREADS; threads += NTHREADS) {
        int argc_2;
        char** argv_2;
        int result;
        snprintf(outputfile, max_path_length, "%s/decode_2%s.sam", TMPDIR, threads ? "threads" : "");
        snprintf(metricsfile, max_path_length, "%s/decode_2%s.metrics", TMPDIR, threads ? "threads" : "");
        snprintf(chksumfile, max_path_length, "%s/decode_2%s.chksum", TMPDIR, threads ? "threads" : "");
        setup_test_2(&argc_2, &argv_2, outputfile, metricsfile, chksumfile, threads);
        main_decode(argc_2-1, argv_2+1);
        free_argv(argc_2,argv_2);

        snprintf(cmd, sizeof(cmd), "diff -I ID:bambi %s %s", outputfile, MKNAME(DATA_DIR,"/out/6383_8_nosplitN.sam"));
        result = system(cmd);
        if (result) {
            fprintf(stderr, "test 2 failed\n");
            failure++;
        } else {
            success++;
        }

        snprintf(cmd, sizeof(cmd), "diff %s %s", chksumfile, MKNAME(DATA_DIR,"/out/decode_2.chksum"));
        result = system(cmd);
        if (result) {
            fprintf(stderr, "test 2 (chksum) failed\n");
            failure++;
        } else {
            success++;
        }
    }

    // check for handling low quality paired reads
    for (int threads = 0; threads <= NTHREADS; threads += NTHREADS) {
        int argc_3;
        char** argv_3;
        int result;
        snprintf(outputfile, max_path_length, "%s/decode_3%s.sam", TMPDIR, threads ? "threads" : "");
        snprintf(metricsfile, max_path_length, "%s/decode_3%s.metrics", TMPDIR, threads ? "threads" : "");
        snprintf(chksumfile, max_path_length, "%s/decode_3%s.chksum", TMPDIR, threads ? "threads" : "");
        setup_test_3(&argc_3, &argv_3, outputfile, metricsfile, chksumfile, threads);
        main_decode(argc_3-1, argv_3+1);
        free_argv(argc_3,argv_3);

        snprintf(cmd, sizeof(cmd), "diff -I ID:bambi %s %s", outputfile, MKNAME(DATA_DIR,"/out/decode_3.sam"));
        result = system(cmd);
        if (result) {
            fprintf(stderr, "test 3 failed\n");
            failure++;
        } else {
            success++;
        }

        snprintf(cmd, sizeof(cmd), "diff %s %s", chksumfile, MKNAME(DATA_DIR,"/out/decode_3.chksum"));
        result = system(cmd);
        if (result) {
            fprintf(stderr, "test 3 (chksum) failed\n");
            failure++;
        } else {
            success++;
        }
    }

    // --dual-tag option
    for (int threads = 0; threads <= NTHREADS; threads += NTHREADS) {
        int argc_4;
        char** argv_4;
        int result;
        snprintf(outputfile, max_path_length,"%s/decode_4%s.sam",TMPDIR, threads ? "threads" : "");
        snprintf(metricsfile, max_path_length, "%s/decode_4%s.metrics", TMPDIR, threads ? "threads" : "");
        setup_test_4(&argc_4, &argv_4, outputfile, metricsfile, threads);
        main_decode(argc_4-1, argv_4+1);
        free_argv(argc_4,argv_4);

        snprintf(cmd, sizeof(cmd), "diff -I ID:bambi %s %s", outputfile, MKNAME(DATA_DIR,"/out/decode_4.sam"));
        result = system(cmd);
        if (result) {
            fprintf(stderr, "test 4 failed at SAM file diff\n");
            failure++;
        } else {
            success++;
        }

        snprintf(cmd, sizeof(cmd), "diff -I ID:bambi %s %s", metricsfile, MKNAME(DATA_DIR,"/out/decode_4.metrics"));
        result = system(cmd);
        if (result) {
            fprintf(stderr, "test 4 failed at metrics file diff\n");
            failure++;
        } else {
            success++;
        }

        snprintf(cmd, sizeof(cmd), "diff -I ID:bambi %s %s", strcat(metricsfile, ".hops"), MKNAME(DATA_DIR,"/out/decode_4.metrics.hops"));
        result = system(cmd);
        if (result) {
            fprintf(stderr, "test 4 failed at tag hops file diff\n");
            failure++;
        } else {
            success++;
        }
    }

    free(metricsfile);
    free(outputfile);
    free(chksumfile);

    printf("decode tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
