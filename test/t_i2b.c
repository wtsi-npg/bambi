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
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>

#include <htslib/kstring.h>

#include "array.h"
#include "bamit.h"

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

ia_t *parseLaneList(char *arg);

int verbose = 0;

int main_i2b(int argc, char *argv[]);

const char * bambi_version(void)
{
    return "12.34";
}

int success = 0;
int failure = 0;

void setup_param_test(int* argc, char*** argv)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160916_miseq_0966_FC/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup("test/data/out/xxx.sam");
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
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160916_miseq_0966_FC/"));
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("teststudy");
    (*argv)[(*argc)++] = strdup("--basecalls-dir");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160916_miseq_0966_FC/Data/Intensities/BaseCalls/"));
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
    (*argv)[(*argc)++] = strdup("16");
    (*argv)[(*argc)++] = strdup("--final-cycle");
    (*argv)[(*argc)++] = strdup("69");
    (*argv)[(*argc)++] = strdup("--final-cycle");
    (*argv)[(*argc)++] = strdup("70");
    (*argv)[(*argc)++] = strdup("--first-index-cycle");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--final-index-cycle");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("-S");

    assert(*argc<100);
}

void setup_simple_test(int* argc, char*** argv, char *outputfile, bool verbose)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160916_miseq_0966_FC/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    (*argv)[(*argc)++] = strdup("Study testStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");

    assert(*argc<100);
}

void setup_multiple_lane_test(int* argc, char*** argv, char *outputfile, bool verbose)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160916_miseq_0966_FC/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("all");
    (*argv)[(*argc)++] = strdup("--first-tile");
    (*argv)[(*argc)++] = strdup("1101");
    (*argv)[(*argc)++] = strdup("--tile-limit");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("Test library");
    (*argv)[(*argc)++] = strdup("--sample-alias");
    (*argv)[(*argc)++] = strdup("Test Sample");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("Study testStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");

    assert(*argc<100);
}

void setup_readgroup_test(int* argc, char*** argv, char *outputfile, bool verbose)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160916_miseq_0966_FC/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
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
    (*argv)[(*argc)++] = strdup("Study testStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");
    (*argv)[(*argc)++] = strdup("--read-group-id");
    (*argv)[(*argc)++] = strdup("6000_1");

    assert(*argc<100);
}

void setup_cyclerange_test(int* argc, char*** argv, char *outputfile, bool verbose)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160919_hiseqx_0557_FC/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--first-tile");
    (*argv)[(*argc)++] = strdup("1101");
    (*argv)[(*argc)++] = strdup("--tile-limit");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("TestLibrary");
    (*argv)[(*argc)++] = strdup("--sample-alias");
    (*argv)[(*argc)++] = strdup("TestSample");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("Study TestStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");
    (*argv)[(*argc)++] = strdup("--first-cycle");
    (*argv)[(*argc)++] = strdup("6");
    (*argv)[(*argc)++] = strdup("--final-cycle");
    (*argv)[(*argc)++] = strdup("10");

    assert(*argc<100);
}

void setup_bcread_test(int* argc, char*** argv, char *outputfile, bool verbose)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160919_nextseq_6230_FC/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--first-tile");
    (*argv)[(*argc)++] = strdup("11101");
    (*argv)[(*argc)++] = strdup("--tile-limit");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("TestLibrary");
    (*argv)[(*argc)++] = strdup("--sample-alias");
    (*argv)[(*argc)++] = strdup("TestSample");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("Study TestStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");
    (*argv)[(*argc)++] = strdup("--first-cycle");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--final-cycle");
    (*argv)[(*argc)++] = strdup("2");
    (*argv)[(*argc)++] = strdup("--first-cycle");
    (*argv)[(*argc)++] = strdup("15");
    (*argv)[(*argc)++] = strdup("--final-cycle");
    (*argv)[(*argc)++] = strdup("19");
    (*argv)[(*argc)++] = strdup("--first-index-cycle");
    (*argv)[(*argc)++] = strdup("30");
    (*argv)[(*argc)++] = strdup("--final-index-cycle");
    (*argv)[(*argc)++] = strdup("31");
    (*argv)[(*argc)++] = strdup("--bc-read");
    (*argv)[(*argc)++] = strdup("2");
    (*argv)[(*argc)++] = strdup("--queue-len");
    (*argv)[(*argc)++] = strdup("200000");

    assert(*argc<100);
}

void setup_dualindex_test(int* argc, char*** argv, char *outputfile, bool verbose)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160919_hiseq4000_7984_FC/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--first-tile");
    (*argv)[(*argc)++] = strdup("1101");
    (*argv)[(*argc)++] = strdup("--tile-limit");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("TestLibrary");
    (*argv)[(*argc)++] = strdup("--sample-alias");
    (*argv)[(*argc)++] = strdup("TestSample");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("Study TestStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");
    (*argv)[(*argc)++] = strdup("--no-filter");
    (*argv)[(*argc)++] = strdup("--barcode-tag");
    (*argv)[(*argc)++] = strdup("tr");
    (*argv)[(*argc)++] = strdup("--quality-tag");
    (*argv)[(*argc)++] = strdup("tq");
    (*argv)[(*argc)++] = strdup("--sec-barcode-tag");
    (*argv)[(*argc)++] = strdup("BC");
    (*argv)[(*argc)++] = strdup("--sec-quality-tag");
    (*argv)[(*argc)++] = strdup("QT");

    assert(*argc<100);
}

void setup_tags_test(int* argc, char*** argv, char *outputfile, bool verbose, bool decode, char *metricsfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160919_hiseq2500_4966_FC/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--first-tile");
    (*argv)[(*argc)++] = strdup("1101");
    (*argv)[(*argc)++] = strdup("--tile-limit");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("TestLibrary");
    (*argv)[(*argc)++] = strdup("--sample-alias");
    (*argv)[(*argc)++] = strdup("TestSample");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("Study TestStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");
    (*argv)[(*argc)++] = strdup("--first-cycle");
    (*argv)[(*argc)++] = strdup("1,30");
    (*argv)[(*argc)++] = strdup("--final-cycle");
    (*argv)[(*argc)++] = strdup("2,32");
    (*argv)[(*argc)++] = strdup("--first-index-cycle");
    (*argv)[(*argc)++] = strdup("3,6,11");
    (*argv)[(*argc)++] = strdup("--final-index-cycle");
    (*argv)[(*argc)++] = strdup("5,10,12");
    (*argv)[(*argc)++] = strdup("--barcode-tag");
    (*argv)[(*argc)++] = strdup("b1,b2,b3");
    (*argv)[(*argc)++] = strdup("--quality-tag");
    (*argv)[(*argc)++] = strdup("q1,q2,q3");
    if (decode) {
        (*argv)[(*argc)++] = strdup("--barcode-file");
        (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160919_hiseq2500_4966_FC/barcodes_i2"));
        (*argv)[(*argc)++] = strdup("--barcode-tag-name");
        (*argv)[(*argc)++] = strdup("b2");
        (*argv)[(*argc)++] = strdup("--metrics-file");
        (*argv)[(*argc)++] = strdup(metricsfile);
    }

    assert(*argc<100);
}

void no_separator_test(int* argc, char*** argv, char *outputfile, bool verbose)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160919_hiseq2500_4966_FC/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--first-tile");
    (*argv)[(*argc)++] = strdup("1101");
    (*argv)[(*argc)++] = strdup("--tile-limit");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("TestLibrary");
    (*argv)[(*argc)++] = strdup("--sample-alias");
    (*argv)[(*argc)++] = strdup("TestSample");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("Study TestStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");
    (*argv)[(*argc)++] = strdup("--first-cycle");
    (*argv)[(*argc)++] = strdup("1,30");
    (*argv)[(*argc)++] = strdup("--final-cycle");
    (*argv)[(*argc)++] = strdup("2,32");
    (*argv)[(*argc)++] = strdup("--first-index-cycle");
    (*argv)[(*argc)++] = strdup("3,6,11");
    (*argv)[(*argc)++] = strdup("--final-index-cycle");
    (*argv)[(*argc)++] = strdup("4,9,12");
    (*argv)[(*argc)++] = strdup("--barcode-tag");
    (*argv)[(*argc)++] = strdup("b1,b2,b1");
    (*argv)[(*argc)++] = strdup("--quality-tag");
    (*argv)[(*argc)++] = strdup("q1,q2,q1");
    (*argv)[(*argc)++] = strdup("-S");

    assert(*argc<100);
}
void separator_test(int* argc, char*** argv, char *outputfile, bool verbose, bool decode, char *metricsfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160919_hiseq2500_4966_FC/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--first-tile");
    (*argv)[(*argc)++] = strdup("1101");
    (*argv)[(*argc)++] = strdup("--tile-limit");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("TestLibrary");
    (*argv)[(*argc)++] = strdup("--sample-alias");
    (*argv)[(*argc)++] = strdup("TestSample");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("Study TestStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");
    (*argv)[(*argc)++] = strdup("--first-cycle");
    (*argv)[(*argc)++] = strdup("1,30");
    (*argv)[(*argc)++] = strdup("--final-cycle");
    (*argv)[(*argc)++] = strdup("2,32");
    (*argv)[(*argc)++] = strdup("--first-index-cycle");
    (*argv)[(*argc)++] = strdup("3,6,11");
    (*argv)[(*argc)++] = strdup("--final-index-cycle");
    (*argv)[(*argc)++] = strdup("4,9,12");
    (*argv)[(*argc)++] = strdup("--barcode-tag");
    (*argv)[(*argc)++] = strdup("b1,b2,b1");
    (*argv)[(*argc)++] = strdup("--quality-tag");
    (*argv)[(*argc)++] = strdup("q1,q2,q1");
    if (decode) {
        (*argv)[(*argc)++] = strdup("--barcode-file");
        (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160919_hiseq2500_4966_FC/barcodes_i1i3_sep"));
        (*argv)[(*argc)++] = strdup("--barcode-tag-name");
        (*argv)[(*argc)++] = strdup("b1");
        (*argv)[(*argc)++] = strdup("--metrics-file");
        (*argv)[(*argc)++] = strdup(metricsfile);
    }

    assert(*argc<100);
}


void consecutive_index_test(int* argc, char*** argv, char *outputfile, bool verbose, bool decode, char *metricsfile)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160919_hiseq2500_4966_FC/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--first-tile");
    (*argv)[(*argc)++] = strdup("1101");
    (*argv)[(*argc)++] = strdup("--tile-limit");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("TestLibrary");
    (*argv)[(*argc)++] = strdup("--sample-alias");
    (*argv)[(*argc)++] = strdup("TestSample");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("Study TestStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");
    (*argv)[(*argc)++] = strdup("--first-cycle");
    (*argv)[(*argc)++] = strdup("1,30");
    (*argv)[(*argc)++] = strdup("--final-cycle");
    (*argv)[(*argc)++] = strdup("2,32");
    (*argv)[(*argc)++] = strdup("--first-index-cycle");
    (*argv)[(*argc)++] = strdup("3,5");
    (*argv)[(*argc)++] = strdup("--final-index-cycle");
    (*argv)[(*argc)++] = strdup("4,7");
    if (decode) {
        (*argv)[(*argc)++] = strdup("--barcode-file");
        (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/160919_hiseq2500_4966_FC/barcodes_ci"));
        (*argv)[(*argc)++] = strdup("--metrics-file");
        (*argv)[(*argc)++] = strdup(metricsfile);
    }

    assert(*argc<100);
}

void novaseq_test(int* argc, char*** argv, char *outputfile, bool verbose)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/novaseq/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--no-filter");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("TestLibrary");
    (*argv)[(*argc)++] = strdup("--sample-alias");
    (*argv)[(*argc)++] = strdup("TestSample");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("Study TestStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");

    assert(*argc<100);
}

void novaseq2_test(int* argc, char*** argv, char *outputfile, bool verbose)
{
    *argc = 0;
    *argv = (char**)calloc(sizeof(char*), 100);
    (*argv)[(*argc)++] = strdup("bambi");
    (*argv)[(*argc)++] = strdup("i2b");
    (*argv)[(*argc)++] = strdup("-i");
    (*argv)[(*argc)++] = strdup(MKNAME(DATA_DIR,"/novaseq_corrupt/Data/Intensities"));
    (*argv)[(*argc)++] = strdup("-o");
    (*argv)[(*argc)++] = strdup(outputfile);
    (*argv)[(*argc)++] = strdup("--lane");
    (*argv)[(*argc)++] = strdup("1");
    (*argv)[(*argc)++] = strdup("--no-filter");
    (*argv)[(*argc)++] = strdup("--library-name");
    (*argv)[(*argc)++] = strdup("TestLibrary");
    (*argv)[(*argc)++] = strdup("--sample-alias");
    (*argv)[(*argc)++] = strdup("TestSample");
    (*argv)[(*argc)++] = strdup("--study-name");
    (*argv)[(*argc)++] = strdup("Study TestStudy");
    (*argv)[(*argc)++] = strdup("--run-start-date");
    (*argv)[(*argc)++] = strdup("2011-03-23T00:00:00+0000");
    (*argv)[(*argc)++] = strdup("--fix-blocks");
    if (verbose) (*argv)[(*argc)++] = strdup("--verbose");

    assert(*argc<100);
}

void free_args(char **argv)
{
    for (int n=0; n<100; n++) free(argv[n]);
    free(argv);
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

void compare_metrics(const char *name, const char *expected, const char *result)
{
    char cmd[1024];
    int res;

    snprintf(cmd, sizeof(cmd), "diff -I ID:bambi '%s' '%s'", expected, result);
    res = system(cmd);
    if (res == 1) {
        fprintf(stderr, "%s : files %s and %s differ\n", name, expected, result);
        failure++;
    } else if (res) {
        fprintf(stderr, "Command \"%s\" failed\n", cmd);
        failure++;
    } else {
        success++;
    }
}

int main(int argc, char**argv)
{
    char template[] = "/tmp/bambi.XXXXXX";
    char *TMPDIR = mkdtemp(template);
    size_t filename_len = strlen(TMPDIR)+64;
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

    if (TMPDIR == NULL) {
        fprintf(stderr,"Can't open temporary directory\n");
        exit(1);
    } else {
        if (verbose) fprintf(stderr,"Created temporary directory: %s\n", TMPDIR);
    }
    char *outputfile = calloc(1, filename_len);
    char *metricsfile = calloc(1, filename_len);

    int argc_1;
    char** argv_1;

    // test parseLaneList()
    ia_t *lanes = parseLaneList("5");
    if (lanes->end != 1) {
        fprintf(stderr,"lanes have %d entries: expected 1\n", lanes->end);
        failure++;
    }
    if (lanes->entries[0] != 5) {
        fprintf(stderr,"lanes[0] is %d: expected 5\n", lanes->entries[0]);
        failure++;
    }
    ia_free(lanes);
    lanes = parseLaneList("1-3,5");
    char *s = ia_join(lanes,",");
    if (strcmp(s, "1,2,3,5")) {
        fprintf(stderr,"Lanes are '%s': expected '1,2,3,5'\n", ia_join(lanes,","));
        failure++;
    }
    free(s);
    ia_free(lanes);

    //
    // simple test
    //

    if (verbose) fprintf(stderr,"\n===> Simple test\n");
    snprintf(outputfile, filename_len, "%s/i2b_1.bam", TMPDIR);
    setup_simple_test(&argc_1, &argv_1, outputfile, verbose);
    main_i2b(argc_1-1, argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/test1.bam"), verbose);
    free_args(argv_1);

    //
    // multiple lane test
    //

    if (verbose) fprintf(stderr,"\n===> Multiple Lane test\n");
    snprintf(outputfile, filename_len, "%s/i2b_m.bam", TMPDIR);
    setup_multiple_lane_test(&argc_1, &argv_1, outputfile, verbose);
    main_i2b(argc_1-1, argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/i2b_m.bam"), verbose);
    free_args(argv_1);

    //
    // Test with non-standard read group ID
    //

    if (verbose) fprintf(stderr,"\n===> Read Group ID test\n");
    snprintf(outputfile, filename_len, "%s/i2b_2.bam", TMPDIR);
    setup_readgroup_test(&argc_1, &argv_1, outputfile, verbose);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/test2.bam"), verbose);
    free_args(argv_1);

    //
    // cycle range test
    //
    if (verbose) fprintf(stderr,"\n===> Cycle Range test\n");
    snprintf(outputfile, filename_len, "%s/i2b_4.bam", TMPDIR);
    setup_cyclerange_test(&argc_1, &argv_1, outputfile, verbose);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/test4.bam"), verbose);
    free_args(argv_1);

    //
    // bc-read test
    //
#if 0
    if (verbose) fprintf(stderr,"\n===> bc-read test\n");
    snprintf(outputfile, filename_len, "%s/i2b_5.bam", TMPDIR);
    setup_bcread_test(&argc_1, &argv_1, outputfile, verbose);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles("BC_READ test", outputfile, MKNAME(DATA_DIR,"/out/test5.bam"));
    free_args(argv_1);
#endif

    //
    // dual index run
    //
    if (verbose) fprintf(stderr,"\n===> Dual Index test\n");
    snprintf(outputfile, filename_len, "%s/i2b_6.bam", TMPDIR);
    setup_dualindex_test(&argc_1, &argv_1, outputfile, verbose);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/test6.bam"), verbose);
    free_args(argv_1);

    //
    // multiple barcode tags run
    //
    if (verbose) fprintf(stderr,"\n===> Multiple Tags test\n");
    snprintf(outputfile, filename_len, "%s/i2b_7.bam", TMPDIR);
    setup_tags_test(&argc_1, &argv_1, outputfile, verbose, false, NULL);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/test7.bam"), verbose);
    free_args(argv_1);

    //
    // multiple barcode tags run with decode
    //
    if (verbose) fprintf(stderr,"\n===> Multiple tags with decode test\n");
    snprintf(outputfile, filename_len, "%s/i2b_7_decode.bam", TMPDIR);
    snprintf(metricsfile, filename_len, "%s/i2b_7_decode.bam.metrics", TMPDIR);
    setup_tags_test(&argc_1, &argv_1, outputfile, verbose, true, metricsfile);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/test7_decode.sam"), verbose);
    compare_metrics("Multiple barcode tags test with decode", MKNAME(DATA_DIR,"/out/test7_decode.bam.metrics"), metricsfile);
    free_args(argv_1);

    //
    // no separator test
    //
    if (verbose) fprintf(stderr,"\n===> no Separator test\n");
    snprintf(outputfile, filename_len, "%s/i2b_8.bam", TMPDIR);
    no_separator_test(&argc_1, &argv_1, outputfile, verbose);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/test8.bam"), verbose);
    free_args(argv_1);

    //
    // separator test
    //
    if (verbose) fprintf(stderr,"\n===> Separator test\n");
    snprintf(outputfile, filename_len, "%s/i2b_9.bam", TMPDIR);
    separator_test(&argc_1, &argv_1, outputfile, verbose, false, NULL);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/test9.bam"), verbose);
    free_args(argv_1);

    //
    // separator test with decode
    //
    if (verbose) fprintf(stderr,"\n===> Separator test with decode\n");
    snprintf(outputfile, filename_len, "%s/i2b_9_decode.bam", TMPDIR);
    snprintf(metricsfile, filename_len, "%s/i2b_9_decode.bam.metrics", TMPDIR);
    separator_test(&argc_1, &argv_1, outputfile, verbose, true, metricsfile);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/test9_decode.sam"), verbose);
    compare_metrics("separator test with decode", MKNAME(DATA_DIR,"/out/test9_decode.bam.metrics"), metricsfile);
    snprintf(metricsfile, filename_len, "%s/i2b_9_decode.bam.metrics.hops", TMPDIR);
    compare_metrics("separator test with decode", MKNAME(DATA_DIR,"/out/test9_decode.bam.metrics.hops"), metricsfile);
    free_args(argv_1);

    //
    // consecutive index test
    //
    if (verbose) fprintf(stderr,"\n===> consecutive test\n");
    snprintf(outputfile, filename_len, "%s/i2b_10.bam", TMPDIR);
    consecutive_index_test(&argc_1, &argv_1, outputfile, verbose, false, NULL);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/test10.bam"), verbose);
    free_args(argv_1);

    //
    // consecutive index test with decode
    //
    if (verbose) fprintf(stderr,"\n===> consecutive test with decode\n");
    snprintf(outputfile, filename_len, "%s/i2b_10.bam", TMPDIR);
    snprintf(metricsfile, filename_len, "%s/i2b_10.bam.metrics", TMPDIR);
    consecutive_index_test(&argc_1, &argv_1, outputfile, verbose, true, metricsfile);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/test10_decode.sam"), verbose);
    compare_metrics("consecutive index test", MKNAME(DATA_DIR,"/out/test10_decode.bam.metrics"), metricsfile);
    snprintf(metricsfile, filename_len, "%s/i2b_10.bam.metrics.hops", TMPDIR);
    compare_metrics("consecutive index test", MKNAME(DATA_DIR,"/out/test10_decode.bam.metrics.hops"), metricsfile);
    free_args(argv_1);

    //
    // novaseq test
    //
    if (verbose) fprintf(stderr,"\n===> NovaSeq test\n");
    snprintf(outputfile, filename_len, "%s/novaseq_1.sam", TMPDIR);
    novaseq_test(&argc_1, &argv_1, outputfile, verbose);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles(outputfile, MKNAME(DATA_DIR,"/out/novaseq_1.sam"), verbose);
    free_args(argv_1);

    //
    // novaseq test with corrupt cbcl file
    //
    if (verbose) fprintf(stderr,"\n===> NovaSeq with corrupt cbcl file test\n");
    snprintf(outputfile, filename_len, "%s/novaseq_2.sam", TMPDIR);
    novaseq2_test(&argc_1, &argv_1, outputfile, verbose);
    main_i2b(argc_1-1,argv_1+1);
    checkFiles("NovaSeq test", outputfile, MKNAME(DATA_DIR,"/out/novaseq_2.sam"));
    free_args(argv_1);

    free(outputfile);
    free(metricsfile);

    printf("i2b tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
