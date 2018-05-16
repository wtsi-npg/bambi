/*  test/t_sf.c -- spatial_filter test cases.

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

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#define xMKNAME(d,f) #d f
#define MKNAME(d,f) xMKNAME(d,f)

int main_spatial_filter(int argc, char *argv[]);

const char * bambi_version(void)
{
    return "12.34";
}

int success = 0;
int failure = 0;

void checkFiles(char *tmpdir, char *gotfile, char *expectfile, int verbose)
{
    char cmd[1024];

    if (verbose) fprintf(stderr,"\nComparing headers: %s with %s\n", gotfile, expectfile);
    // compare headers
    sprintf(cmd,"samtools view -H %s |grep -v ^@PG |sort | perl -n -e 'chomp; @x=split /\t/;@y=sort @x; print join \",\",@y; print \"\n\";' | sed s:/tmp/bambi[^/]*:/tmp/xyzzy:g > %s/got.txt", gotfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"samtools view -H %s | grep -v ^@PG| sort | perl -n -e 'chomp; @x=split /\t/;@y=sort @x; print join \",\",@y; print \"\n\";' | sed s:/tmp/bambi[^/]*:/tmp/xyzzy:g > %s/expect.txt", expectfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"diff %s/got.txt %s/expect.txt", tmpdir, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }

    // compare records
    if (verbose) fprintf(stderr,"\nComparing records: %s with %s\n", gotfile, expectfile);
    sprintf(cmd,"samtools view %s > %s/got.txt", gotfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"samtools view %s > %s/expect.txt", expectfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"diff %s/got.txt %s/expect.txt", tmpdir, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
}

void checkFilterFiles(char *prog, char *tmpdir, char *gotfile, char *expectfile)
{
    char cmd[1024];
    sprintf(cmd,"%s -D -v -F %s > %s/got.txt", prog, gotfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"%s -D -v -F %s > %s/expect.txt", prog, expectfile, tmpdir);
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    sprintf(cmd,"diff -ICommand %s/got.txt %s/expect.txt", tmpdir, tmpdir);
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
    char outputfile[512];
    char filterfile[512];
    char cmd[2048];
    char prog[512];

    sprintf(prog, "%s", "src/bambi spatial_filter");

    // create filter
    if (verbose) fprintf(stderr,"Creating filter\n");
    snprintf(filterfile, sizeof(filterfile), "%s/sf_1.filter", TMPDIR);
    snprintf(outputfile, sizeof(outputfile), "%s/sf_filtered.bam", TMPDIR);
    snprintf(cmd, sizeof(cmd), "%s -c -F %s %s", prog, filterfile, MKNAME(DATA_DIR,"/sf.bam"));
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    checkFilterFiles(prog, TMPDIR, filterfile, MKNAME(DATA_DIR,"/out/sf_1.filter"));

    // apply filter
    if (verbose) fprintf(stderr,"Applying filter\n");
    snprintf(cmd, sizeof(cmd), "%s -a --verbose -F %s -o %s %s", prog, filterfile, outputfile, MKNAME(DATA_DIR,"/sf.bam"));
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    checkFiles(TMPDIR, outputfile, MKNAME(DATA_DIR,"/out/sf_filtered.bam"), verbose);

    // multiple filters
    if (verbose) fprintf(stderr,"Multiple filters\n");
    snprintf(outputfile, sizeof(outputfile), "%s/sf_filtered.bam", TMPDIR);
    snprintf(cmd, sizeof(cmd), "%s -a -v -f -F %s,%s,%s --rg 25077_3#3,25077_4#3,25077_5#3 -o %s %s", prog, MKNAME(DATA_DIR,"/sf_1.filter"), MKNAME(DATA_DIR,"/sf_2.filter"), MKNAME(DATA_DIR,"/sf_3.filter"), outputfile, MKNAME(DATA_DIR,"/sf2.bam"));
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    checkFiles(TMPDIR, outputfile, MKNAME(DATA_DIR,"/out/sf2.bam"), verbose);

    // missing filters
    if (verbose) fprintf(stderr,"Missing filters\n");
    snprintf(outputfile, sizeof(outputfile), "%s/sf_filtered_2.bam", TMPDIR);
    snprintf(cmd, sizeof(cmd), "%s -a -v -f -F %s,%s,%s --rg 25077_3#3,25077_4#3,25077_5#x -o %s %s", prog, MKNAME(DATA_DIR,"/sf_1.filter"), MKNAME(DATA_DIR,"/sf_2.filter"), MKNAME(DATA_DIR,"/sf_3.filter"), outputfile, MKNAME(DATA_DIR,"/sf2.bam"));
    if (system(cmd)) { fprintf(stderr,"Command failed: %s\n",cmd); failure++; }
    checkFiles(TMPDIR, outputfile, MKNAME(DATA_DIR,"/out/sf_filtered_2.bam"), verbose);

    printf("spatial_filter tests: %s\n", failure ? "FAILED" : "Passed");
    return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
