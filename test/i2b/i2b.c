/*  test/i2b/i2b.c -- i2b test cases.

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
#include <config.h>

#include "bambi.h"
#include "../../i2b.c"
#include <stdlib.h>
#include <unistd.h>
#include "version.h"

int verbose = 0;

const char * bambi_version(void)
{
    return "12.34";
}

int success = 0;
int failure = 0;

void setup_test_1(int* argc, char*** argv)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup("test/i2b/110323_HS13_06000_B_B039WABXX/Data/Intensities");
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup("test/i2b/out/xxx.sam");
    (*argv)[(*argc)++] = strdup("--output-fmt");
    (*argv)[(*argc)++] = strdup("sam");
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--first-tile");
    (*argv)[(*argc)++] = strdup("1101");
    (*argv)[(*argc)++] = strdup("--tile-limit");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("testlibrary");
    (*argv)[(*argc)++] = strdup("--run-folder");
    (*argv)[(*argc)++] = strdup("test/i2b/110323_HS13_06000_B_B039WABXX/");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("teststudy");
    (*argv)[(*argc)++] = strdup("--basecalls-dir");
    (*argv)[(*argc)++] = strdup("test/i2b/110323_HS13_06000_B_B039WABXX/Data/Intensities/BaseCalls/");
    (*argv)[(*argc)++] = strdup("--generate-secondary-basecalls");
    (*argv)[(*argc)++] = strdup("--no-filter");
    (*argv)[(*argc)++] = strdup("--sequencing-centre");
    (*argv)[(*argc)++] = strdup("XY");
    (*argv)[(*argc)++] = strdup("--platform");
    (*argv)[(*argc)++] = strdup("Illumina");
    (*argv)[(*argc)++] = strdup("--first-tile");
    (*argv)[(*argc)++] = strdup("1103");
    (*argv)[(*argc)++] = strdup("--tile-limit");
    (*argv)[(*argc)++] = strdup("5");
    (*argv)[(*argc)++] = strdup("--barcode-tag");
    (*argv)[(*argc)++] = strdup("AB");
    (*argv)[(*argc)++] = strdup("--quality-tag");
    (*argv)[(*argc)++] = strdup("CD");
    (*argv)[(*argc)++] = strdup("--sec-barcode-tag");
    (*argv)[(*argc)++] = strdup("WX");
    (*argv)[(*argc)++] = strdup("--sec-quality-tag");
    (*argv)[(*argc)++] = strdup("YZ");
    (*argv)[(*argc)++] = strdup("--bc-read");
    (*argv)[(*argc)++] = strdup("2");
    (*argv)[(*argc)++] = strdup("--first-cycle");
    (*argv)[(*argc)++] = strdup("7");
    (*argv)[(*argc)++] = strdup("--first-cycle");
    (*argv)[(*argc)++] = strdup("17");
    (*argv)[(*argc)++] = strdup("--first-cycle");
    (*argv)[(*argc)++] = strdup("70");
    (*argv)[(*argc)++] = strdup("--final-cycle");
    (*argv)[(*argc)++] = strdup("9");
    (*argv)[(*argc)++] = strdup("--final-index-cycle");
    (*argv)[(*argc)++] = strdup("9");
    (*argv)[(*argc)++] = strdup("--add-cluster-index-tag");

    assert(*argc<100);
}

void setup_test_2(int* argc, char*** argv)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup("test/i2b/110323_HS13_06000_B_B039WABXX/Data/Intensities");
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup("test/i2b/out/xxx.bam");
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--first-tile");
    (*argv)[(*argc)++] = strdup("1101");
    (*argv)[(*argc)++] = strdup("--tile-limit");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("Test library");
    (*argv)[(*argc)++] = strdup("--sample-alias");
    (*argv)[(*argc)++] = strdup("Test Sample");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("testStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    (*argv)[(*argc)++] = strdup("--verbose");

    assert(*argc<100);
}

void checkLike(char *name, char *expected, char *actual)
{
    if (actual == NULL) actual = "<null>";
    if (strstr(actual, expected) == NULL) {
        fprintf(stderr, "%s\n" "Expected: %s\n" "Got:      %s\n", name, expected, actual);
        failure++;
    }
}

void checkEqual(char *name, char *expected, char *actual)
{
    if (actual == NULL) actual = "<null>";
    if (strcmp(expected, actual)) {
        fprintf(stderr, "%s\n" "Expected: %s\n" "Got:      %s\n", name, expected, actual);
        failure++;
    }
}

void icheckEqual(char *name, int expected, int actual)
{
    if (expected != actual) {
        fprintf(stderr, "%s\n" "Expected: %d\n" "Got:      %d\n", name, expected, actual);
        failure++;
    }
}

void test_paramaters(void)
{
    // minimal options
    int argc_1;
    char** argv_1;
    setup_test_1(&argc_1, &argv_1);
    opts_t *opts = i2b_parse_args(argc_1-1, argv_1+1);

    if (verbose) printf("Testing paramaters\n");

    if (!opts) {
        fprintf(stderr, "parse_args failed\n");
        failure++;
        return;
    }

    checkLike("options: intensity-dir", "test/i2b/110323_HS13_06000_B_B039WABXX/Data/Intensities", opts->intensity_dir);
    checkEqual("options: output-file", "test/i2b/out/xxx.sam", opts->output_file);
    checkEqual("options: sample-alias", "testlibrary", opts->sample_alias);
    checkEqual("options: study-name", "teststudy", opts->study_name);
    checkLike("options: run-folder", "test/i2b/110323_HS13_06000_B_B039WABXX", opts->run_folder);
    checkLike("options: basecalls-dir", "test/i2b/110323_HS13_06000_B_B039WABXX/Data/Intensities/BaseCalls", opts->basecalls_dir);
    icheckEqual("options: lane", 1, opts->lane);
    icheckEqual("options: generate-secondary-basecalls", 1, opts->generate_secondary_basecalls);
    icheckEqual("options: no-filter", 1, opts->no_filter);
    checkEqual("options: read-group-id", "1", opts->read_group_id);
    checkEqual("options: sequencing-centre", "XY", opts->sequencing_centre);
    checkEqual("options: platform", "Illumina", opts->platform);
    icheckEqual("options: first-tile", 1103, opts->first_tile);
    icheckEqual("options: tile-limit", 5, opts->tile_limit);
    checkEqual("options: barcode-tag", "AB", opts->barcode_tag);
    checkEqual("options: quality-tag", "CD", opts->quality_tag);
    checkEqual("options: sec-barcode-tag", "WX", opts->barcode_tag2);
    checkEqual("options: sec-quality-tag", "YZ", opts->quality_tag2);
    icheckEqual("options: sec-bc-read", 2, opts->sec_bc_read);
    icheckEqual("options: first-cycle", 3, opts->first_cycle->end);
    icheckEqual("options: first-cycle[0]", 7, opts->first_cycle->entries[0]);
    icheckEqual("options: first-cycle[1]", 17, opts->first_cycle->entries[1]);
    icheckEqual("options: first-cycle[2]", 70, opts->first_cycle->entries[2]);
    icheckEqual("options: final-cycle", 1, opts->final_cycle->end);
    icheckEqual("options: final-index-cycle", 1, opts->final_index_cycle->end);
    icheckEqual("options: final-cycle[0]", 9, opts->final_cycle->entries[0]);
    icheckEqual("options: add-cluster-index-tag", 1, opts->add_cluster_index_tag);
}

int main(int argc, char**argv)
{
    int result;

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

    test_paramaters();


    int argc_2;
    char** argv_2;
    setup_test_2(&argc_2, &argv_2);
    main_i2b(argc_2-1, argv_2+1);

    result = system("diff test/i2b/out/xxx.bam test/i2b/out/test1.bam");
    if (result) {
        fprintf(stderr, "test 2 failed\n");
        failure++;
    } else {
        success++;
    }

    printf("i2b tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
