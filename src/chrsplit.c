/*  chrsplit.c -- split a BAM file into separate files by chromosome.

    Copyright (C) 2016 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "bambi.h"
#include <assert.h>
#include <ctype.h>
#include <htslib/sam.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <unistd.h>
#include <regex.h>
#include <libgen.h>
#include <time.h>
#include <fcntl.h>

#include <cram/sam_header.h>

#include "array.h"
#include "bamit.h"

char *strptime(const char *s, const char *format, struct tm *tm);

static void freeRecord(void *r) { bam_destroy1((bam1_t *)r); }

/*
 * structure to hold options
 */
typedef struct {
    int verbose;
    char *argv_list;
    va_t *subset;
    char compression_level;
    char *in_file;
    char *target_file;
    char *exclude_file;
    char *output_fmt;
    char *input_fmt;
    bool exclude_unaligned;
    bool invert;
} opts_t;

/*
 * Release all the options
 */

static void free_opts(opts_t* opts)
{
    if (!opts) return;
    free(opts->argv_list);
    free(opts->in_file);
    free(opts->target_file);
    free(opts->exclude_file);
    free(opts->output_fmt);
    free(opts->input_fmt);
    va_free(opts->subset);
    free(opts);
}

/*
 * display usage information
 */
static void usage(FILE *write_to)
{
    fprintf(write_to,
"Usage: bambi chrsplit [options]\n"
"\n"
"Program to split a BAM file based on given chromosomes.  More precisely:\n"
"Define subsets of the @SQ reference sequences in the BAM header by their SN fields.\n"
"Filter pairs of reads based on alignment records referring to a member of this @SQ subset.\n"
"Specify chromosome subset on the command line or use the default of MT and Y\n"
"Send pairs of reads for which either read aligns to a chromosome/@SQ not specified in the given set to an 'excluded' output.\n"
"Send other pairs, i.e. both unaligned, both aligned to the given subset, or one aligned to the given subset and the other\n"
"unaligned to a 'target'\n"
"\n"
"Options:\n"
"  -i   --input                 BAM file to read\n"
"  -o   --output                BAM file to for target reads\n"
"  -e   --exclude               BAM file to for excluded reads\n"
"  -s   --subset                Reference sequences to target. [Default 'Y,MT']\n"
"  -u   --exclude-unaligned     Exclude read groups where all reads are unaligned\n"
"  -V   --invert                Treat the -s option as a list to exclude rather than target\n"
"  -v   --verbose               verbose output\n"
"       --input-fmt             [sam/bam/cram] [default: bam]\n"
"       --output-fmt            [sam/bam/cram] [default: bam]\n"
"       --compression-level     [0..9]\n"
);
}

/*
 * Takes the command line options and turns them into something we can understand
 */
opts_t* chrsplit_parse_args(int argc, char *argv[])
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char* optstring = "vi:o:n:e:s:Vu";

    static const struct option lopts[] = {
        { "verbose",            0, 0, 'v' },
        { "input",              1, 0, 'i' },
        { "output",             1, 0, 'o' },
        { "exclude",            1, 0, 'e' },
        { "subset",             1, 0, 's' },
        { "exclude-unaligned",  0, 0, 'u' },
        { "invert",             0, 0, 'V' },
        { "compression-level",  1, 0, 0 },
        { "input-fmt",          1, 0, 0 },
        { "output-fmt",         1, 0, 0 },
        { NULL, 0, NULL, 0 }
    };

    opts_t* opts = calloc(sizeof(opts_t), 1);
    if (!opts) { perror("cannot allocate option parsing memory"); return NULL; }

    opts->argv_list = stringify_argv(argc+1, argv-1);
    if (opts->argv_list[strlen(opts->argv_list)-1] == ' ') opts->argv_list[strlen(opts->argv_list)-1] = 0;

    opts->subset = va_init(5,free);

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, optstring, lopts, &option_index)) != -1) {
        const char *arg;
        switch (opt) {
        case 's':   parse_tags(opts->subset, optarg);
                    break;
        case 'i':   opts->in_file = strdup(optarg);
                    break;
        case 'o':   opts->target_file = strdup(optarg);
                    break;
        case 'e':   opts->exclude_file = strdup(optarg);
                    break;
        case 'v':   opts->verbose++;
                    break;
        case 'u':   opts->exclude_unaligned = true;
                    break;
        case 'V':   opts->invert = true;
                    break;
        case 0:     arg = lopts[option_index].name;
                         if (strcmp(arg, "output-fmt") == 0)              opts->output_fmt = strdup(optarg);
                    else if (strcmp(arg, "input-fmt") == 0)               opts->input_fmt = strdup(optarg);
                    else if (strcmp(arg, "compression-level") == 0)       opts->compression_level = *optarg;
                    else {
                        fprintf(stderr,"\nUnknown option: %s\n\n", arg); 
                        usage(stdout); free_opts(opts);
                        return NULL;
                    }
                    break;
        default:    fprintf(stderr,"Unknown option: '%c'\n", opt);
            /* else fall-through */
        case '?':   usage(stdout); free_opts(opts); return NULL;
        }
    }

    argc -= optind;
    argv += optind;

    optind = 0;

    if (!opts->in_file) {
        fprintf(stderr,"You must specify an input file\n");
        usage(stderr); return NULL;
    }

    if (!opts->target_file) {
        fprintf(stderr,"You must specify a target file\n");
        usage(stderr); return NULL;
    }

    if (!opts->exclude_file) {
        fprintf(stderr,"You must specify an exclude file\n");
        usage(stderr); return NULL;
    }

    if (opts->compression_level && !isdigit(opts->compression_level)) {
        fprintf(stderr, "compression-level must be a digit in the range [0..9], not '%c'\n", opts->compression_level);
        usage(stderr); return NULL;
    }

    if (va_isEmpty(opts->subset)) {
        va_push(opts->subset,strdup("MT"));
        va_push(opts->subset,strdup("Y"));
    }

    return opts;
}

char *getReferenceName(bam1_t *rec, bam_hdr_t *h)
{
    if (rec->core.tid == -1) return NULL;
    return h->target_name[rec->core.tid];
}

static int find_in_subset(va_t *va, char *search)
{
    int n;
    for (n=0; n < va->end; n++) {
        if (strcmp((char *)va->entries[n],search) == 0) return n;
    }
    return -1;
}

/*
 * convert SAM_hdr to bam_hdr
 */
static void sam_hdr_unparse(SAM_hdr *sh, bam_hdr_t *h)
{
    free(h->text);
    sam_hdr_rebuild(sh);
    h->text = strdup(sam_hdr_str(sh));
    h->l_text = sam_hdr_length(sh);
    sam_hdr_free(sh);
}

/*
 * open a single BAM file
 */
static samFile *openSamFile(char *fname, char *fmt, char compression, char rw)
{
    samFile *f = NULL;
    htsFormat *format = NULL;
    char mode[] = "xbC";

    if (fmt) {
        format = calloc(1,sizeof(htsFormat));
        if (hts_parse_format(format, fmt) < 0) {
            fprintf(stderr,"Unknown input format: %s\n", fmt);
            exit(1);
        }
    }
    mode[0] = rw;
    mode[2] = compression;
    f = hts_open_format(fname, mode, format);
    free(format);
    if (!f) {
        fprintf(stderr,"Could not open file (%s)\n", fname);
        exit(1);
    }

    return f;
}

/*
 * add PG line to output BAM file
 */
static void addPGLine(BAMit_t *bit, opts_t *opts, char *ot)
{
    SAM_hdr *sh = sam_hdr_parse_(bit->h->text,bit->h->l_text);

    // specify sort order
    sh->sort_order = ORDER_UNSORTED;

    // add new PG line
    sam_hdr_add_PG(sh, "bambi",
                   "OT", ot,
                   "VN", bambi_version(),
                   "CL", opts->argv_list,
                   "DS", "Split BAM file by chromosomes",
                   NULL, NULL);
    sam_hdr_unparse(sh,bit->h);
}

/*
 * Read records from a given iterator until the qname changes
 */
static va_t *read_record_set(BAMit_t *bit, char *qname)
{
    va_t *recordSet = va_init(5,freeRecord);

    while (BAMit_hasnext(bit) && strcmp(bam_get_qname(BAMit_peek(bit)),qname) == 0) {
        bam1_t *rec = bam_init1();
        bam_copy1(rec,BAMit_next(bit));
        va_push(recordSet,rec);
    }

    return recordSet;
}

/*
 * write a record set to an output BAM file
 */
static void writeRecordSet(BAMit_t *bit, va_t *recordSet)
{
    int n,r;

    for (n=0; n < recordSet->end; n++) {
        bam1_t *rec = recordSet->entries[n];
        r = sam_write1(bit->f, bit->h, rec);
        if (r <= 0) {
            fprintf(stderr, "Problem writing record %d : %d\n", n,r);
            exit(1);
        }
    }
}

/*
 * Read the input file, write to the output files
 */
static int processFiles(BAMit_t *in_bam, BAMit_t *target_bam, BAMit_t *exclude_bam, opts_t *opts)
{
    BAMit_t *outBam;
    int n;
    bool all_unaligned;

    addPGLine(target_bam,opts,"TARGET");
    addPGLine(exclude_bam,opts,"EXCLUDED");
    if (sam_hdr_write(target_bam->f, target_bam->h)) { fprintf(stderr,"Failed to write target header\n"); exit(1); }
    if (sam_hdr_write(exclude_bam->f, target_bam->h)) { fprintf(stderr,"Failed to write exclude header\n"); exit(1); }

    while (BAMit_hasnext(in_bam)) {
        // for each record set
        bam1_t *rec = BAMit_peek(in_bam);
        char *qname = strdup(bam_get_qname(rec));
        va_t *recordSet = read_record_set(in_bam,qname);
        free(qname);

        outBam = target_bam;
        all_unaligned = true;

        for (n=0; n < recordSet->end; n++) {
            bam1_t *rec = recordSet->entries[n];
            if (!(rec->core.flag & BAM_FUNMAP)) {
                bool notfound = (find_in_subset(opts->subset,getReferenceName(rec,in_bam->h))==-1);
                if (( opts->invert) ^ (find_in_subset(opts->subset,getReferenceName(rec,in_bam->h))==-1))
                {
                    outBam = exclude_bam;
                    break;
                } else {
                    all_unaligned = false;
                }
            } else { // one or both reads are unmapped
                if (opts->exclude_unaligned) {
                    outBam = exclude_bam;
                }
            }
        }
        if (all_unaligned && opts->exclude_unaligned) outBam = exclude_bam;

        writeRecordSet(outBam,recordSet);
        va_free(recordSet);
    }

    return 0;
}

/*
 * Main code
 *
 * Open and read a BAM file, and split into two output files
 */
static int chrsplit(opts_t* opts)
{
    int retcode = 1;

    samFile *f = NULL;
    bam_hdr_t *h = NULL;

    BAMit_t *in_bam = NULL;
    BAMit_t *target_bam = NULL;
    BAMit_t *exclude_bam = NULL;

    f = openSamFile(opts->in_file, opts->input_fmt, opts->compression_level, 'r');
    h = sam_hdr_read(f);
    in_bam = BAMit_init(f,h);

    f = openSamFile(opts->target_file, opts->output_fmt, opts->compression_level, 'w');
    target_bam = BAMit_init(f,bam_hdr_dup(h));

    f = openSamFile(opts->exclude_file, opts->output_fmt, opts->compression_level, 'w');
    exclude_bam = BAMit_init(f,bam_hdr_dup(h));

    if (in_bam && target_bam && exclude_bam) {
        retcode = processFiles(in_bam, target_bam, exclude_bam, opts);
    }

    // tidy up after us
    BAMit_free(in_bam);
    BAMit_free(target_bam);
    BAMit_free(exclude_bam);

    return retcode;
}

/*
 * called from bambi to perform selection by alignment
 *
 * Parse the command line arguments, then call the main select() function
 *
 * returns 0 on success, 1 if there was a problem
 */
int main_chrsplit(int argc, char *argv[])
{
    int ret = 1;
    opts_t* opts = chrsplit_parse_args(argc, argv);
    if (opts) ret = chrsplit(opts);
    free_opts(opts);
    return ret;
}
