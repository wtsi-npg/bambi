/*  bambi.c -- main bambi command front-end.

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

#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#include "htslib/hts.h"
#include "bambi.h"

int main_decode(int argc, char *argv[]);
int main_i2b(int argc, char *argv[]);
int main_select(int argc, char *argv[]);

const char *bambi_version()
{
    return VERSION;
}

static void vprint_error_core(const char *subcommand, const char *format, va_list args, const char *extra)
{
    fflush(stdout);
    if (subcommand && *subcommand) fprintf(stderr, "bambi %s: ", subcommand);
    else fprintf(stderr, "bambi: ");
    vfprintf(stderr, format, args);
    if (extra) fprintf(stderr, ": %s\n", extra);
    else fprintf(stderr, "\n");
    fflush(stderr);
}

void print_error(const char *subcommand, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    vprint_error_core(subcommand, format, args, NULL);
    va_end(args);
}

void print_error_errno(const char *subcommand, const char *format, ...)
{
    int err = errno;
    va_list args;
    va_start(args, format);
    vprint_error_core(subcommand, format, args, strerror(err));
    va_end(args);
}

static void usage(FILE *fp)
{
    /* Please improve the grouping */

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
    else if (strcmp(argv[1], "--version") == 0) {
        printf( "bambi %s\n"
                "Using htslib %s\n"
                "Copyright (C) 2016 Genome Research Ltd.\n",
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
