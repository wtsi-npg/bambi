/*  test/decode/decode.c -- decode test cases.

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
#include "../src/hts_addendum.c"
#include "../src/decode.c"

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

const char * bambi_version(void)
{
    return "12.34";
}

int success = 0;
int failure = 0;

void setup_test_1(int* argc, char*** argv, char *outputfile)
{
    *argc = 16;
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
    (*argv)[13] = strdup(MKNAME(DATA_DIR,"/out/decode_1.metrics"));
    (*argv)[14] = strdup("--barcode-tag-name");
    (*argv)[15] = strdup("RT");
}

void setup_test_2(int* argc, char*** argv, char *outputfile)
{
    *argc = 18;
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
    (*argv)[15] = strdup(MKNAME(DATA_DIR,"/out/decode_2.metrics"));
    (*argv)[16] = strdup("--barcode-tag-name");
    (*argv)[17] = strdup("RT");
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
    if ((n=countMismatches(a,b)) == e) success++;
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
    test_countMismatches("ABCNXYZ","ABCxXYZ",0);
    test_countMismatches("ABCiXYZ","ABCNXYZ",0);
    test_countMismatches("NBCiXYZ",".BCNXYz",1);
    test_countMismatches("AGCACGTT","AxCACGTTXXXXXX",1);

    //
    // Now test the actual decoding
    //

    // minimal options
    int argc_1;
    char** argv_1;
    char *outputfile = calloc(1,strlen(TMPDIR)+64);
    sprintf(outputfile,"%s/decode_1.sam", TMPDIR);

    setup_test_1(&argc_1, &argv_1, outputfile);
    main_decode(argc_1-1, argv_1+1);

    char *cmd = calloc(1,1024);
    sprintf(cmd,"diff -I ID:bambi %s %s", outputfile, MKNAME(DATA_DIR,"/out/6383_9_nosplit_nochange.sam"));
    int result = system(cmd);
    if (result) {
        fprintf(stderr, "test 1 failed\n");
        failure++;
    } else {
        success++;
    }

    // --convert_low_quality option
    int argc_2;
    char** argv_2;
    sprintf(outputfile,"%s/decode_2.sam",TMPDIR);
    setup_test_2(&argc_2, &argv_2, outputfile);
    main_decode(argc_2-1, argv_2+1);

    sprintf(cmd,"diff -I ID:bambi %s %s", outputfile, MKNAME(DATA_DIR,"/out/6383_8_nosplitN.sam"));
    result = system(cmd);
    if (result) {
        fprintf(stderr, "test 2 failed\n");
        failure++;
    } else {
        success++;
    }

    printf("decode tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
