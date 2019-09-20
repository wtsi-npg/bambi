/*  bambi.c -- main bambi command front-end.

    Copyright (C) 2016-2019 Genome Research Ltd.

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

#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#include "htslib/hts.h"
#include "bambi.h"
#include "bambi_utils.h"

int main_decode(int argc, char *argv[]);
int main_i2b(int argc, char *argv[]);
int main_select(int argc, char *argv[]);
int main_chrsplit(int argc, char *argv[]);
int main_read2tags(int argc, char *argv[]);
int main_spatial_filter(int argc, char *argv[]);
int main_seqchksum(int argc, char *argv[]);
int main_adapters(int argc, char *argv[]);
int main_update(int argc, char *argv[]);

const char *bambi_version()
{
    return VERSION;
}

static void usage(FILE *fp)
{
    fprintf(fp,
"\n"
"Program: bambi (Tools for alignments in the SAM format)\n"
"Version: %s (using htslib %s)\n\n", bambi_version(), hts_version());
    fprintf(fp,
"Usage:   bambi <command> [options]\n"
"\n"
"Commands:\n"
"     decode         decode a multiplexed SAM/BAM/CRAM file by read groups\n"
"     i2b            converts illumina files to SAM/BAM/CRAM files\n"
"     select         select reads by alignment\n"
"     chrsplit       split reads by chromosome\n"
"     read2tags      convert reads into tags\n"
"     spatial_filter spatial filtering\n"
"     seqchksum      calculate checksums for a bam file\n"
"     adapters       find and remove adapters\n"
"     update         update an existing BAM/SAM/CRAM file\n"
"\n"
"bambi <command> for help on a particular command\n"
"\n");
}

int main(int argc, char *argv[])
{
    if (argc < 2) { usage(stderr); return 1; }

    if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0) {
        if (argc == 2) { usage(stdout); return 0; }

        // Otherwise change "bambi help COMMAND [...]" to "bambi COMMAND";
        // main_xyz() functions by convention display the subcommand's usage
        // when invoked without any arguments.
        argv++;
        argc = 2;
    }

    int ret = 0;
         if (strcmp(argv[1], "decode") == 0)    ret = main_decode(argc-1, argv+1);
    else if (strcmp(argv[1], "i2b") == 0)       ret = main_i2b(argc-1, argv+1);
    else if (strcmp(argv[1], "select") == 0)    ret = main_select(argc-1, argv+1);
    else if (strcmp(argv[1], "chrsplit") == 0)  ret = main_chrsplit(argc-1, argv+1);
    else if (strcmp(argv[1], "read2tags") == 0) ret = main_read2tags(argc-1, argv+1);
    else if (strcmp(argv[1], "spatial_filter") == 0) ret = main_spatial_filter(argc-1, argv+1);
    else if (strcmp(argv[1], "seqchksum") == 0) ret = main_seqchksum(argc-1, argv+1);
    else if (strcmp(argv[1], "adapters") == 0) ret = main_adapters(argc-1, argv+1);
    else if (strcmp(argv[1], "update") == 0) ret = main_update(argc-1, argv+1);
    else if (strcmp(argv[1], "--version") == 0) {
        printf( "bambi %s\n"
                "Using htslib %s\n"
                "Copyright (C) 2017 Genome Research Ltd.\n",
                bambi_version(), hts_version());
    }
    else if (strcmp(argv[1], "--version-only") == 0) {
        printf("%s+htslib-%s\n", bambi_version(), hts_version());
    }
    else {
        fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
        usage(stderr);
        return 1;
    }
    return ret;
}
