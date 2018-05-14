/*  decode.c -- index decoder subcommand.

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

#include <bambi.h>
#include <assert.h>
#include <htslib/sam.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <unistd.h>
#include <regex.h>
#include <errno.h>
#include <htslib/khash.h>
#include <htslib/thread_pool.h>
#include <cram/sam_header.h>
#include <inttypes.h>

#include "decode.h"
#include "bamit.h"
#include "hash_table.h"

#define xstr(s) str(s)
#define str(s) #s

#define DEFAULT_MAX_LOW_QUALITY_TO_CONVERT 15
#define DEFAULT_MAX_NO_CALLS 2
#define DEFAULT_MAX_MISMATCHES 1
#define DEFAULT_MIN_MISMATCH_DELTA 1
#define DEFAULT_BARCODE_TAG "BC"
#define DEFAULT_QUALITY_TAG "QT"
#define TEMPLATES_PER_JOB 5000

// Size of stack allocations to use for storing barcodes.  If too small, malloc will be used instead.
// Ideally this should be bigger than the longest barcode expected.
#ifndef STACK_BC_LEN
#define STACK_BC_LEN 32
#endif

enum match {
    MATCHED_NONE,
    MATCHED_FIRST,
    MATCHED_SECOND,
    MATCHED_BOTH,
    MATCHED_NEW
};


/*
 * structure to hold options
 */
struct decode_opts_t {
    char *input_name;
    char *output_name;
    char *barcode_name;
    char *metrics_name;
    char *barcode_tag_name;
    char *quality_tag_name;
    bool verbose;
    int max_low_quality_to_convert;
    bool convert_low_quality;
    int max_no_calls;
    int max_mismatches;
    int min_mismatch_delta;
    bool change_read_name;
    char *argv_list;
    char *input_fmt;
    char *output_fmt;
    char compression_level;
    int nthreads;
    int idx1_len, idx2_len;
    bool ignore_pf;
    unsigned short dual_tag;
};

decode_opts_t *decode_init_opts(int argc, char **argv)
{
    decode_opts_t *opts = calloc(1, sizeof(*opts));
    if (!opts) die("Out of memory");
    opts->argv_list = stringify_argv(argc, argv);
    if (opts->argv_list[strlen(opts->argv_list)-1] == ' ') opts->argv_list[strlen(opts->argv_list)-1] = 0;

    // set defaults
    opts->max_low_quality_to_convert = DEFAULT_MAX_LOW_QUALITY_TO_CONVERT;
    opts->max_no_calls = DEFAULT_MAX_NO_CALLS;
    opts->max_mismatches = DEFAULT_MAX_MISMATCHES;
    opts->min_mismatch_delta = DEFAULT_MIN_MISMATCH_DELTA;
    opts->verbose = false;
    opts->convert_low_quality = false;
    opts->change_read_name = false;
    opts->barcode_tag_name = NULL;
    opts->quality_tag_name = NULL;
    opts->ignore_pf = 0;
    opts->dual_tag = 0;
    return opts;
}

void decode_free_opts(decode_opts_t *opts)
{
    if (!opts) return;
    free(opts->input_name);
    free(opts->output_name);
    free(opts->barcode_name);
    free(opts->barcode_tag_name);
    free(opts->quality_tag_name);
    free(opts->argv_list);
    free(opts->input_fmt);
    free(opts->output_fmt);
    free(opts->metrics_name);
    free(opts);
}

/*
 * details read from barcode file
 * Plus metrics information for each barcode
 */
typedef struct {
    char *seq, *idx1, *idx2;
    char *name;
    char *lib;
    char *sample;
    char *desc;
    uint64_t reads, pf_reads, perfect, pf_perfect, one_mismatch, pf_one_mismatch;
} bc_details_t;

// Data for thread pool jobs
typedef struct decode_thread_data_t {
    va_t *record_set;                   // records to process
    ia_t *template_counts;              // records in each template
    va_t *barcode_array;                // job-local copy of barcodes array
    bc_details_t *barcodes;             // memory for job-local barcodes
    HashTable *tagHopHash;              // job-local tag hops hash
    HashTable *barcodeHash;             // pointer to shared barcodeHash
    decode_opts_t *opts;                       // pointer to shared opts
    int nrec;                           // number of live records in record_set
    int result;                         // job result, 0 = success
    struct decode_thread_data_t *next;  // for free list
} decode_thread_data_t;

static bc_details_t *bcd_init(void)
{
    bc_details_t *bcd = calloc(1, sizeof(bc_details_t));
    bcd->seq = NULL;
    bcd->idx1 = NULL;
    bcd->idx2 = NULL;
    bcd->name = NULL;
    bcd->lib = NULL;
    bcd->sample = NULL;
    bcd->desc = NULL;
    return bcd;
}

/*
 * Print metrics file header
 */
static void print_header(FILE* f, decode_opts_t* opts, bool metrics) {
    // print header
    fprintf(f, "##\n");
    fprintf(f, "# BARCODE_TAG_NAME=%s ", opts->barcode_tag_name);
    fprintf(f, "MAX_MISMATCHES=%d ", opts->max_mismatches);
    fprintf(f, "MIN_MISMATCH_DELTA=%d ", opts->min_mismatch_delta);
    fprintf(f, "MAX_NO_CALLS=%d ", opts->max_no_calls);
    fprintf(f, "\n");
    fprintf(f, "##\n");
    fprintf(f, "# ID:bambi VN:%s (htslib %s) CL:%s\n", bambi_version(), hts_version(), opts->argv_list);
    fprintf(f, "\n");
    fprintf(f, "##\n");
    fprintf(f, "BARCODE\t");
    if (metrics) {
        fprintf(f, "BARCODE_NAME\t");
        fprintf(f, "LIBRARY_NAME\t");
        fprintf(f, "SAMPLE_NAME\t");
        fprintf(f, "DESCRIPTION\t");
    }
    fprintf(f, "READS\t");
    if (!opts->ignore_pf) {
        fprintf(f, "PF_READS\t");
    }
    fprintf(f, "PERFECT_MATCHES\t");
    if (!opts->ignore_pf) {
        fprintf(f, "PF_PERFECT_MATCHES\t");
    }
    if (metrics) {
        fprintf(f, "ONE_MISMATCH_MATCHES\t");
        if (!opts->ignore_pf) {
            fprintf(f, "PF_ONE_MISMATCH_MATCHES\t");
        }
    }
    fprintf(f, "PCT_MATCHES\t");
    fprintf(f, "RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT");
    if (!opts->ignore_pf) {
        fprintf(f, "\tPF_PCT_MATCHES");
    }
    if (!opts->ignore_pf) {
        fprintf(f, "\tPF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT");
    }
    if (!opts->ignore_pf) {
        fprintf(f, "\tPF_NORMALIZED_MATCHES");
    }
    fprintf(f, "\n");
}

void free_bcd(void *entry)
{
    bc_details_t *bcd = (bc_details_t *)entry;
    free(bcd->seq);
    free(bcd->idx1);
    free(bcd->idx2);
    free(bcd->name);
    free(bcd->lib);
    free(bcd->sample);
    free(bcd->desc);
    free(bcd);
}

void free_taghop_bcd(void *entry)
{
    bc_details_t *bcd = (bc_details_t *)entry;
    free(bcd->seq);
    free(bcd);
}

static int compareTagHops(const void *t1, const void *t2) {
    bc_details_t *th1 = *(bc_details_t **)t1;
    bc_details_t *th2 = *(bc_details_t **)t2;

    int read_diff = th1->reads - th2->reads;
    //if read count is equal, sort by number of perfect matches
    return read_diff ? -read_diff : -(th1->perfect - th2->perfect);
}

static void sortTagHops(va_t *tagHopArray) {
    qsort(tagHopArray->entries, tagHopArray->end, sizeof(bc_details_t*), compareTagHops);
}

static void freeRecord(void *r) { bam_destroy1((bam1_t *)r); }

/*
 * display usage information
 */
static void usage(FILE *write_to)
{
    fprintf(write_to,
"Usage: bambi decode [options] filename\n"
"\n"
"Options:\n"
"  -o   --output                        output file [default: stdout]\n"
"  -v   --verbose                       verbose output\n"
"  -b   --barcode-file                  file containing barcodes\n"
"       --convert-low-quality           Convert low quality bases in barcode read to 'N'\n"
"       --max-low-quality-to-convert    Max low quality phred value to convert bases in barcode\n"
"                                       read to 'N' [default: " xstr(DEFAULT_MAX_LOW_QUALITY_TO_CONVERT) "]\n"
"       --max-no-calls                  Max allowable number of no-calls in a barcode read before\n"
"                                       it is considered unmatchable [default: " xstr(DEFAULT_MAX_NO_CALLS) "]\n"
"       --max-mismatches                Maximum mismatches for a barcode to be considered a match\n"
"                                       [default: " xstr(DEFAULT_MAX_MISMATCHES) "]\n"
"       --min-mismatch-delta            Minimum difference between number of mismatches in the best\n"
"                                       and second best barcodes for a barcode to be considered a\n"
"                                       match [default: " xstr(DEFAULT_MIN_MISMATCH_DELTA) "]\n"
"       --change-read-name              Change the read name by adding #<barcode> suffix\n"
"       --metrics-file                  Per-barcode and per-lane metrics written to this file\n"
"       --barcode-tag-name              Barcode tag name [default: " DEFAULT_BARCODE_TAG "]\n"
"       --quality-tag-name              Quality tag name [default: " DEFAULT_QUALITY_TAG "]\n"
"       --input-fmt                     format of input file [sam/bam/cram]\n"
"       --output-fmt                    format of output file [sam/bam/cram]\n"
"       --compression-level             Compression level of output file [0..9]\n"
"  -t   --threads                       number of threads to use [default: 1]\n"
"       --ignore-pf                     Doesn't output PF statistics\n"
"       --dual-tag                      Dual tag position in the barcode string (between 2 and barcode length - 1)\n"
);
}

/*
 * Functions to set some options (used by i2b)
 */
void set_decode_opt_barcode_name(decode_opts_t *opts, const char *name)
{
    free(opts->barcode_name);
    opts->barcode_name = strdup(name);
    if (!opts->barcode_name) die("Out of memory");
}

void set_decode_opt_metrics_name(decode_opts_t *opts, const char *name)
{
    free(opts->metrics_name);
    opts->metrics_name = strdup(name);
    if (!opts->metrics_name) die("Out of memory");
}

void set_decode_opt_barcode_tag_name(decode_opts_t *opts, const char *name) {
    free(opts->barcode_tag_name);
    opts->barcode_tag_name = strdup(name);
    if (!opts->barcode_tag_name) die("Out of memory");
}

void set_decode_opt_max_low_quality_to_convert(decode_opts_t *opts, int val)
{
    opts->max_low_quality_to_convert = val;
}

void set_decode_opt_convert_low_quality(decode_opts_t *opts, bool flag)
{
    opts->convert_low_quality = flag;
}

void set_decode_opt_max_no_calls(decode_opts_t *opts, int val)
{
    opts->max_no_calls = val;
}

void set_decode_opt_max_mismatches(decode_opts_t *opts, int val)
{
    opts->max_mismatches = val;
}

void set_decode_opt_min_mismatch_delta(decode_opts_t *opts, int val)
{
    opts->min_mismatch_delta = val;
}

void set_decode_opt_change_read_name(decode_opts_t *opts, bool flag)
{
    opts->change_read_name = flag;
}

void set_decode_opt_ignore_pf(decode_opts_t *opts, bool flag)
{
    opts->ignore_pf = flag;
}

/*
 * Takes the command line options and turns them into something we can understand
 */
static decode_opts_t* parse_args(int argc, char *argv[])
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char* optstring = "i:o:vb:t:";

    static const struct option lopts[] = {
        { "input",                      1, 0, 'i' },
        { "output",                     1, 0, 'o' },
        { "verbose",                    0, 0, 'v' },
        { "max-low-quality-to-convert", 1, 0, 0 },
        { "convert-low-quality",        0, 0, 0 },
        { "barcode-file",               1, 0, 'b' },
        { "max-no-calls",               1, 0, 0 },
        { "max-mismatches",             1, 0, 0 },
        { "min-mismatch-delta",         1, 0, 0 },
        { "change-read-name",           0, 0, 0 },
        { "metrics-file",               1, 0, 0 },
        { "barcode-tag-name",           1, 0, 0 },
        { "quality-tag-name",           1, 0, 0 },
        { "input-fmt",                  1, 0, 0 },
        { "output-fmt",                 1, 0, 0 },
        { "compression-level",          1, 0, 0 },
        { "ignore-pf",                  0, 0, 0 },
        { "dual-tag",                   1, 0, 0 },
        { "threads",                    1, 0, 't' },
        { NULL, 0, NULL, 0 }
    };

    decode_opts_t* opts = decode_init_opts(argc+1, argv-1);
    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, optstring, lopts, &option_index)) != -1) {
        const char *arg;
        switch (opt) {
        case 'i':   opts->input_name = strdup(optarg);
                    break;
        case 'o':   opts->output_name = strdup(optarg);
                    break;
        case 'v':   opts->verbose = true;
                    break;
        case 'b':   opts->barcode_name = strdup(optarg);
                    break;
        case 't':   opts->nthreads = atoi(optarg);
                    break;
        case 0:     arg = lopts[option_index].name;
                         if (strcmp(arg, "metrics-file") == 0)               opts->metrics_name = strdup(optarg);
                    else if (strcmp(arg, "max-low-quality-to-convert") == 0) opts->max_low_quality_to_convert = atoi(optarg);
                    else if (strcmp(arg, "convert-low-quality") == 0)        opts->convert_low_quality = true;
                    else if (strcmp(arg, "max-no-calls") == 0)               opts->max_no_calls = atoi(optarg);
                    else if (strcmp(arg, "max-mismatches") == 0)             opts->max_mismatches = atoi(optarg);
                    else if (strcmp(arg, "min-mismatch-delta") == 0)         opts->min_mismatch_delta = atoi(optarg);
                    else if (strcmp(arg, "change-read-name") == 0)           opts->change_read_name = true;
                    else if (strcmp(arg, "barcode-tag-name") == 0)           opts->barcode_tag_name = strdup(optarg);
                    else if (strcmp(arg, "quality-tag-name") == 0)           opts->quality_tag_name = strdup(optarg);
                    else if (strcmp(arg, "input-fmt") == 0)                  opts->input_fmt = strdup(optarg);
                    else if (strcmp(arg, "output-fmt") == 0)                 opts->output_fmt = strdup(optarg);
                    else if (strcmp(arg, "compression-level") == 0)          opts->compression_level = *optarg;
                    else if (strcmp(arg, "ignore-pf") == 0)                  opts->ignore_pf = true;
                    else if (strcmp(arg, "dual-tag") == 0)                  {opts->dual_tag = (short)atoi(optarg);
                                                                             opts->max_no_calls = 0;}  
                    else {
                        printf("\nUnknown option: %s\n\n", arg); 
                        usage(stdout); decode_free_opts(opts);
                        return NULL;
                    }
                    break;
        default:    printf("Unknown option: '%c'\n", opt);
            /* else fall-through */
        case '?':   usage(stdout); decode_free_opts(opts); return NULL;
        }
    }

    argc -= optind;
    argv += optind;

    if (argc > 0) opts->input_name = strdup(argv[0]);
    optind = 0;

    // some validation and tidying
    if (!opts->input_name) {
        fprintf(stderr,"You must specify an input file (-i or --input)\n");
        usage(stderr); decode_free_opts(opts);
        return NULL;
    }
    if (!opts->barcode_name) {
        fprintf(stderr,"You must specify a barcode (tags) file (-b or --barcode-file)\n");
        usage(stderr); decode_free_opts(opts);
        return NULL;
    }

    if (!opts->barcode_tag_name) opts->barcode_tag_name = strdup(DEFAULT_BARCODE_TAG);
    if (!opts->quality_tag_name) opts->quality_tag_name = strdup(DEFAULT_QUALITY_TAG);

    // output defaults to stdout
    if (!opts->output_name) opts->output_name = strdup("-");

    return opts;
}

//
// checkBarcodeQuality(char *newBarcode, char *barcode, size_t barcode_len, char *quality, size_t quality_len, decode_opts_t *opts);
//
// Write a new barcode read string with low quality bases converted to 'N' into newBarcode.
//
static void checkBarcodeQuality(char *newBarcode,
                                char *bc_tag, size_t bc_len,
                                char *qt_tag, size_t qt_len, decode_opts_t *opts)
{
    memcpy(newBarcode, bc_tag, bc_len + 1);
    if (!qt_tag) return;

    if (bc_len != qt_len) {
        die("checkBarcodeQuality(): barcode and quality are different lengths\n");
    }

    int mlq = opts->max_low_quality_to_convert ? opts->max_low_quality_to_convert 
                                               : DEFAULT_MAX_LOW_QUALITY_TO_CONVERT;
    for (size_t i=0; i < qt_len; i++) {
        int qual = qt_tag[i] - 33;
        if (isalpha(newBarcode[i]) && (qual <= mlq)) newBarcode[i] = 'N';
    }
}

void writeMetricsLine(FILE *f, bc_details_t *bcd, decode_opts_t *opts, uint64_t total_reads, uint64_t max_reads, uint64_t total_pf_reads, uint64_t max_pf_reads, uint64_t total_pf_reads_assigned, uint64_t nReads, bool metrics)
{
    fprintf(f, "%s", bcd->idx1);
    if (bcd->idx2 && *bcd->idx2) fprintf(f, "-%s", bcd->idx2);
    fprintf(f, "\t");
    if (metrics) {
        fprintf(f, "%s\t", bcd->name);
        fprintf(f, "%s\t", bcd->lib);
        fprintf(f, "%s\t", bcd->sample);
        fprintf(f, "%s\t", bcd->desc);
    }
    fprintf(f, "%"PRIu64"\t", bcd->reads);
    if (!opts->ignore_pf) {
        fprintf(f, "%"PRIu64"\t", bcd->pf_reads); 
    }
    fprintf(f, "%"PRIu64"\t", bcd->perfect);
    if (!opts->ignore_pf) {
        fprintf(f, "%"PRIu64"\t", bcd->pf_perfect); 
    }
    if (metrics) {
        fprintf(f, "%"PRIu64"\t", bcd->one_mismatch);
        if (!opts->ignore_pf) {
            fprintf(f, "%"PRIu64"\t", bcd->pf_one_mismatch); 
        }
    }
    fprintf(f, "%.3f\t", total_reads ? bcd->reads / (double)total_reads : 0 );
    fprintf(f, "%.3f", max_reads ? bcd->reads / (double)max_reads  : 0 );
    if (!opts->ignore_pf) {
        fprintf(f, "\t%.3f", total_pf_reads ? bcd->pf_reads / (double)total_pf_reads  : 0 ); 
    }
    if (!opts->ignore_pf) {
        fprintf(f, "\t%.3f", max_pf_reads ? bcd->pf_reads / (double)max_pf_reads  : 0 ); 
    }
    if (!opts->ignore_pf) {
        fprintf(f, "\t%.3f", total_pf_reads_assigned ? bcd->pf_reads * nReads / (double)total_pf_reads_assigned  : 0);
    }
    fprintf(f, "\n");
}


/*
 *
 */
int writeMetrics(va_t *barcodeArray, HashTable *tagHopHash, decode_opts_t *opts)
{
    bc_details_t *bcd = barcodeArray->entries[0];
    uint64_t total_reads = bcd->reads;
    uint64_t total_pf_reads = bcd->pf_reads;
    uint64_t total_pf_reads_assigned = 0;
    uint64_t max_reads = 0;
    uint64_t total_original_reads = 0;
    uint64_t total_hop_reads = 0;
    uint64_t max_pf_reads = 0;
    uint64_t nReads = 0;
    int n;

    // Open the metrics file
    FILE *f = fopen(opts->metrics_name, "w");
    if (!f) {
        fprintf(stderr,"Can't open metrics file %s\n", opts->metrics_name);
        return 1;
    }

    // first loop to count things
    for (n=1; n < barcodeArray->end; n++) {
        bc_details_t *bcd = barcodeArray->entries[n];;
        total_reads += bcd->reads;
        total_original_reads += bcd->reads;
        total_pf_reads += bcd->pf_reads;
        total_pf_reads_assigned += bcd->pf_reads;
        if (max_reads < bcd->reads) max_reads = bcd->reads;
        if (max_pf_reads < bcd->pf_reads) max_pf_reads = bcd->pf_reads;
        nReads++;
    }

    // Copy the tag hop hash into an array and sort it
    va_t *tagHopArray = NULL;
    if (tagHopHash && tagHopHash->nused) {
        tagHopArray = va_init(barcodeArray->end, NULL);
        HashIter *iter = HashTableIterCreate();
        HashItem *hi;
        while ( (hi = HashTableIterNext(tagHopHash, iter)) != NULL) {
            va_push(tagHopArray, hi->data.p);
        }
        HashTableIterDestroy(iter);
        sortTagHops(tagHopArray);
    }

    if (tagHopArray) {
        for (n=0; n < tagHopArray->end; n++) {
            bc_details_t *bcd = tagHopArray->entries[n];
            total_hop_reads += bcd->reads;
        }
    }

    // print header
    print_header(f, opts, true);

    // second loop to print things
    for (n=1; n < barcodeArray->end; n++) {
        bc_details_t *bcd = barcodeArray->entries[n];
        writeMetricsLine(f, bcd, opts, total_reads, max_reads, total_pf_reads, max_pf_reads, total_pf_reads_assigned, nReads, true);
    }
    // treat Tag 0 as a special case
    bcd = barcodeArray->entries[0];
    bcd->perfect = 0;
    bcd->pf_perfect = 0;
    bcd->name[0] = 0;
    writeMetricsLine(f, bcd, opts, total_reads, max_reads, total_pf_reads, max_pf_reads, 0, nReads, true);

    fclose(f);

    /*
     * Now write tag hop metrics file - if there are any
     */

    if (opts->idx2_len) {
        char *metrics_hops_name = malloc(strlen(opts->metrics_name)+6);
        strcpy(metrics_hops_name, opts->metrics_name);
        strcat(metrics_hops_name, ".hops");
        FILE *g = fopen(metrics_hops_name, "w");
        if (!g) {
            fprintf(stderr,"Can't open tag hops file %s\n", metrics_hops_name);
        } else {
            fprintf(g, "##\n");
            fprintf(g, "# TOTAL_READS=%"PRIu64", ", total_reads);
            fprintf(g, "TOTAL_ORIGINAL_TAG_READS=%"PRIu64", ", total_original_reads);
            fprintf(g, "TOTAL_TAG_HOP_READS=%"PRIu64", ", total_hop_reads);
            fprintf(g, "MAX_READ_ON_A_TAG=%"PRIu64", ", max_reads);
            fprintf(g, "TOTAL_TAG_HOPS=%d, ", (tagHopArray ? tagHopArray->end : 0));
            fprintf(g, "PCT_TAG_HOPS=%f\n", total_reads ? ((float)total_hop_reads / total_reads * 100) : 0.00) ;
            print_header(g, opts, false);

            if (tagHopArray) {
                sortTagHops(tagHopArray);
                for (n=0; n < tagHopArray->end; n++) {
                    bc_details_t *bcd = tagHopArray->entries[n];
                    writeMetricsLine(g, bcd, opts, total_reads, max_reads, total_pf_reads, max_pf_reads, total_pf_reads_assigned, nReads, false);
                }
            }
            fclose(g);
        }
        free(metrics_hops_name);
    }

    va_free(tagHopArray);
    return 0;
}

/*
 * split a dual index (eg ACACAC-TGTGTG) into two different indexes.
 * If a single index is given, then the second index is an empty string.
 */
static void split_index(char *seq, size_t seq_len, int dual_tag, char **idx1_ptr, char **idx2_ptr, size_t idx1_sz, size_t idx2_sz)
{
    size_t idx1_len, idx2_start, idx2_len;
    if (dual_tag) {
        idx1_len = dual_tag - 1;
        idx2_start = dual_tag;
        idx2_len = seq_len - dual_tag;
    } else {
        idx1_len = strcspn(seq, INDEX_SEPARATOR);
        idx2_start = idx1_len + strspn(seq + idx1_len, INDEX_SEPARATOR);
        idx2_len = strcspn(seq + idx2_start, INDEX_SEPARATOR);
    }
    if (idx1_len >= idx1_sz) {
        *idx1_ptr = malloc(idx1_len + 1);
        if (!*idx1_ptr) die("Out of memory");
    }
    memcpy(*idx1_ptr, seq, idx1_len); (*idx1_ptr)[idx1_len] = '\0';
    if (idx2_len >= idx2_sz) {
        *idx2_ptr = malloc(idx2_len + 1);
        if (!*idx2_ptr) die("Out of memory");
    }
    memcpy(*idx2_ptr, seq + idx2_start, idx2_len); (*idx2_ptr)[idx2_len] = '\0';
}

/*
 * Read the barcode file into an array
 */
va_t *loadBarcodeFile(decode_opts_t *opts)
{
    int lineno = 0;
    int idx1_len=0, idx2_len=0;
    va_t *barcodeArray = va_init(100,free_bcd);

    // initialise first entry for null metrics
    bc_details_t *bcd = bcd_init();
    bcd->name   = strdup("0");
    bcd->lib    = strdup("");
    bcd->sample = strdup("");
    bcd->desc   = strdup("");
    va_push(barcodeArray,bcd);

    FILE *fh = fopen(opts->barcode_name,"r");
    if (!fh) {
        fprintf(stderr,"ERROR: Can't open barcode file %s\n", opts->barcode_name);
        return NULL;
    }
    
    char *buf = NULL;
    size_t n;
    lineno++;
    if (getline(&buf,&n,fh) < 0) {    // burn first line which is a header
        fprintf(stderr,"ERROR: problem reading barcode file\n");
        return NULL;
    }
    free(buf); buf=NULL;

    while (getline(&buf, &n, fh) > 0) {
        lineno++;
        char *s;
        if (buf[strlen(buf)-1] == '\n') buf[strlen(buf)-1]=0;   // remove trailing lf
        bc_details_t *bcd = bcd_init();
        s = strtok(buf,"\t");  if (!s) die("Can't read sequence from tag file: Line %d\n", lineno);
        bcd->seq     = strdup(s);
        s = strtok(NULL,"\t"); if (!s) die("Can't read name from tag file: Line %d\n", lineno);
        bcd->name    = strdup(s);
        s = strtok(NULL,"\t"); bcd->lib     = s ? strdup(s) : strdup("");
        s = strtok(NULL,"\t"); bcd->sample  = s ? strdup(s) : strdup("");
        s = strtok(NULL,"\t"); bcd->desc    = s ? strdup(s) : strdup("");

        split_index(bcd->seq, strlen(bcd->seq), opts->dual_tag, &bcd->idx1, &bcd->idx2, 0, 0);

        va_push(barcodeArray,bcd);
        free(buf); buf=NULL;

        if (idx1_len == 0) {
            idx1_len = strlen(bcd->idx1);
            idx2_len = strlen(bcd->idx2);
        } else {
            if ( (idx1_len != strlen(bcd->idx1)) && (idx2_len != strlen(bcd->idx2)) ) {
                fprintf(stderr,"ERROR: Tag '%s' is a different length to the previous tag\n", bcd->seq);
                return NULL;
            }
        }
    }

    opts->idx1_len = idx1_len;
    opts->idx2_len = idx2_len;
    bcd = barcodeArray->entries[0];
    bcd->idx1 = calloc(1,idx1_len+1); memset(bcd->idx1, 'N', idx1_len);
    bcd->idx2 = calloc(1,idx2_len+1); memset(bcd->idx2, 'N', idx2_len);
    bcd->seq = calloc(1,idx1_len+idx2_len+2);
    strcpy(bcd->seq,bcd->idx1);
    if (idx2_len) strcat(bcd->seq,INDEX_SEPARATOR);
    strcat(bcd->seq,bcd->idx2);

    free(buf);
    fclose(fh);
    return barcodeArray;
}

/*
 * return true if base is a noCall
 */
//#define isNoCall(c) (c=='N')
int isNoCall(char b)
{
    return b=='N' || b=='n' || b=='.';
}

/*
 * Count the number of noCalls in a sequence
 */
static int noCalls(char *s)
{
    int n=0;
    while (*s) {
        if (isNoCall(*s++)) n++;
    }
    return n;
}

/*
 * count number of mismatches between two sequences
 * (ignoring noCalls)
 */
static int countMismatches(char *tag, char *barcode, int maxval)
{
    int n = 0;
    for (int i=0; tag[i] && barcode[i]; i++) {
        if ((tag[i] != barcode[i]) && (barcode[i] != 'N')) {
            n++;
            if (n>maxval) return n;     // exit early if we can
        }
    }
    return n;
}

/*
 * For a failed match, check is there is tag hopping to report
 */
static bc_details_t *check_tag_hopping(char *barcode, va_t *barcodeArray, HashTable *tagHopHash, decode_opts_t *opts)
{
    bc_details_t *bcd = NULL, *best_match1 = NULL, *best_match2 = NULL;
    char stack_idx1[STACK_BC_LEN], stack_idx2[STACK_BC_LEN];
    char *idx1 = stack_idx1, *idx2 = stack_idx2;
    int nmBest1 = opts->idx1_len + opts->idx2_len + 1;
    int nmBest2 = nmBest1;

    split_index(barcode, strlen(barcode), opts->dual_tag, &idx1, &idx2, sizeof(stack_idx1), sizeof(stack_idx2));

    // for each tag in barcodeArray
    for (int n=1; n < barcodeArray->end; n++) {
        bc_details_t *bcd = barcodeArray->entries[n];

        int nMismatches1 = countMismatches(bcd->idx1, idx1, nmBest1);
        int nMismatches2 = countMismatches(bcd->idx2, idx2, nmBest2);

        // match the first tag
        if (nMismatches1 < nmBest1) {
            nmBest1 = nMismatches1;
            best_match1 = bcd;
        }

        // match the second tag
        if (nMismatches2 < nmBest2) {
            nmBest2 = nMismatches2;
            best_match2 = bcd;
        }
    }

    if (idx1 != stack_idx1) free(idx1);
    if (idx2 != stack_idx2) free(idx2);

    bool matched_first = (nmBest1 == 0 );
    bool matched_second = (nmBest2 == 0 );

    if (matched_first && matched_second) {
        HashData hd;
        HashItem *hi;
        char stack_key[STACK_BC_LEN];
        char *key = stack_key;

        assert(best_match1 != NULL && best_match2 != NULL);
        if (opts->idx1_len + opts->idx2_len + 2 >= sizeof(stack_key)) {
            key = malloc(opts->idx1_len + opts->idx2_len + 2);
            if (!key) die("Out of memory");
        }
        memcpy(key, best_match1->idx1, opts->idx1_len);
        memcpy(key + opts->idx1_len, INDEX_SEPARATOR, 1);
        memcpy(key + opts->idx1_len + 1, best_match2->idx2, opts->idx2_len);
        key[opts->idx1_len + opts->idx2_len + 1] = '\0';
        hi = HashTableSearch(tagHopHash, key, 0);
        if (hi) {
            bcd = hi->data.p;
        } else {
            bcd = calloc(1, sizeof(bc_details_t)); //create a new entry with the two tags
            bcd->idx1 = best_match1->idx1;
            bcd->idx2 = best_match2->idx2;
            bcd->seq = strdup(key);
            bcd->name = "0";
            bcd->lib = "DUMMY_LIB";
            bcd->sample = "DUMMY_SAMPLE";
            bcd->desc = NULL;
            hd.p = bcd;
            HashTableAdd(tagHopHash, key, 0, hd, NULL);
        }
        if (key != stack_key) free(key);
    }

    return bcd;
}


/*
 * find the best match in the barcode (tag) file for a given barcode
 * return the tag, if a match found, else return NULL
 */
bc_details_t *findBestMatch(char *barcode, va_t *barcodeArray, HashTable *barcodeHash, decode_opts_t *opts)
{
    int bcLen = opts->idx1_len + opts->idx2_len + 1;   // size of barcode sequence in barcode file
    bc_details_t *best_match = NULL;
    int nmBest = bcLen;             // number of mismatches (best)
    int nm2Best = bcLen;            // number of mismatches (second best)

    bool matched = false;

    // First look in the barcodeHash for an exact match
    // This optimisation only applies when mismatch_delta is 1 (the default)
    if (opts->min_mismatch_delta <= 1) {
        HashItem *hi;
        hi = HashTableSearch(barcodeHash, barcode, 0);
        if (hi) {
            return barcodeArray->entries[hi->data.i];
        }
    }

    // No exact match, so do it the hard way...
    for (int n=1; n < barcodeArray->end; n++) {
        bc_details_t *bcd = barcodeArray->entries[n];

        int nMismatches = countMismatches(bcd->seq, barcode, nm2Best);
        if (nMismatches < nmBest) {
            nm2Best = nmBest;
            nmBest = nMismatches;
            best_match = bcd;
        } else {
            if (nMismatches < nm2Best) nm2Best = nMismatches;
        }
    }

    matched = best_match && nmBest <= opts->max_mismatches && nm2Best - nmBest >= opts->min_mismatch_delta;

    if (!matched) {
        best_match = barcodeArray->entries[0];
    }

    return best_match;
}

/*
 * Update the metrics information
 */
static void updateMetrics(bc_details_t *bcd, char *seq, bool isPf)
{
    int n = 99;
    if (seq) n = countMismatches(bcd->seq, seq, 999);

    bcd->reads++;
    if (isPf) bcd->pf_reads++;

    if (n==0) {     // count perfect matches
        bcd->perfect++;
        if (isPf) bcd->pf_perfect++;
    }

    if (n==1) {     // count out-by-one matches
        bcd->one_mismatch++;
        if (isPf) bcd->pf_one_mismatch++;
    }
}

/*
 * find the best match in the barcode (tag) file, and return the corresponding barcode name
 * If no match found, check for tag hopping, and return dummy entry 0
 */
char *findBarcodeName(char *barcode, va_t *barcodeArray, HashTable *barcodeHash, HashTable *tagHopHash, decode_opts_t *opts, bool isPf, bool isUpdateMetrics)
{
    bc_details_t *bcd;
    if (noCalls(barcode) > opts->max_no_calls) {
        bcd = barcodeArray->entries[0];
        if (isUpdateMetrics) updateMetrics(bcd, barcode, isPf);
    } else {
        bcd = findBestMatch(barcode, barcodeArray, barcodeHash, opts);
        if (isUpdateMetrics) updateMetrics(bcd, barcode, isPf);
        if ((bcd == barcodeArray->entries[0]) && opts->idx2_len) {
            bc_details_t *tag_hop = check_tag_hopping(barcode, barcodeArray, tagHopHash, opts);
            if (isUpdateMetrics && tag_hop) updateMetrics(tag_hop, barcode, isPf);
        }
    }
    return bcd->name;
}

/*
 * make a new tag by appending #<name> to the old tag
 */
static void makeNewTag(bam1_t *rec, char *tag, char *name, char **newtag, size_t newtag_sz)
{
    char *rg = "";
    uint8_t *p = bam_aux_get(rec,tag);
    size_t name_len = strlen(name), rg_len;
    if (p) rg = bam_aux2Z(p);
    rg_len = strlen(rg);
    if (name_len + rg_len + 2 > newtag_sz) {
        *newtag = malloc(name_len + rg_len + 2);
        if (!*newtag) die("Out of memory");
    }
    memcpy(*newtag, rg, rg_len);
    (*newtag)[rg_len] = '#';
    memcpy(*newtag + rg_len + 1, name, name_len + 1);
}

/*
 * Change the read name by adding "#<suffix>"
 */
static void add_suffix(bam1_t *rec, char *suffix)
{
    int oldqlen = strlen((char *)rec->data);
    int newlen = rec->l_data + strlen(suffix) + 1;

    if (newlen > rec->m_data) {
        rec->m_data = newlen;
        kroundup32(rec->m_data);
        rec->data = (uint8_t *)realloc(rec->data, rec->m_data);
    }
    memmove(rec->data + oldqlen + strlen(suffix) + 1,
            rec->data + oldqlen,
            rec->l_data - oldqlen);
    rec->data[oldqlen] = '#';
    memmove(rec->data + oldqlen + 1,
            suffix,
            strlen(suffix) + 1);
    rec->l_data = newlen;
    rec->core.l_qname += strlen(suffix) + 1;
}

/*
 * Add a new @RG line to the header
 */
static void addNewRG(SAM_hdr *sh, char *entry, char *bcname, char *lib, char *sample, char *desc)
{
    assert(entry != NULL);
    char *saveptr;
    char *p = strtok_r(entry,"\t",&saveptr);
    char *newtag = malloc(strlen(p)+1+strlen(bcname)+1);
    strcpy(newtag, p);
    strcat(newtag,"#");
    strcat(newtag, bcname);
    sam_hdr_add(sh, "RG", "ID", newtag, NULL, NULL);

    SAM_hdr_type *hdr = sam_hdr_find(sh, "RG", "ID", newtag);
    while (1) {
        char *pu = NULL;
        char *t = strtok_r(NULL, ":", &saveptr);
        if (!t) break;
        char *v = strtok_r(NULL, "\t", &saveptr);
        if (!v) break;

        // handle special cases
        if (strcmp(t,"PU") == 0) {
            // add #bcname
            pu = malloc(strlen(v) + 1 + strlen(bcname) + 1);
            strcpy(pu, v); strcat(pu,"#"); strcat(pu,bcname);
            v = pu;
        }
        if (strcmp(t,"LB") == 0) {
            if (lib) v = lib;        // use library name
        }
        if (strcmp(t,"DS") == 0) {
            if (desc) v = desc;       // use desc
        }
        if (strcmp(t,"SM") == 0) {
            if (sample) v = sample;     // use sample name
        }
        sam_hdr_update(sh, hdr, t, v, NULL);
        if (pu) free(pu);
    }
    free(newtag);
}

/*
 * for each "@RG ID:x" in the header, replace with
 * "@RG IDx#barcode" for each barcode
 *
 * And don't forget to add a @PG header
 */ 
static void changeHeader(va_t *barcodeArray, bam_hdr_t *h, char *argv_list)
{
    SAM_hdr *sh = sam_hdr_parse_(h->text, h->l_text);
    char **rgArray = malloc(sizeof(char*) * sh->nrg);
    int nrg = sh->nrg;
    int i, n;

    sam_hdr_add_PG(sh, "bambi", "VN", bambi_version(), "CL", argv_list, NULL);

    // store the RG names
    for (n=0; n < sh->nrg; n++) {
        // store the names and tags as a string <name>:<tag>:<val>:<tag>:<val>...
        // eg 1:PL:Illumina:PU:110608_HS19

        // first pass to determine size of string required
        int sz=strlen(sh->rg[n].name)+1;
        SAM_hdr_tag *rgtag = sh->rg[n].tag;
        while (rgtag) {
            if (strncmp(rgtag->str,"ID:",3)) {  // ignore name
                sz += 3 + strlen(rgtag->str) + 1;
            }
            rgtag = rgtag->next;
        }
        char *entry = malloc(sz+1);

        // second pass to create string
        strcpy(entry,sh->rg[n].name);
        rgtag = sh->rg[n].tag;
        while (rgtag) {
            if (strncmp(rgtag->str,"ID:",3)) {  // ignore name
                strcat(entry,"\t");
                strcat(entry,rgtag->str);
            }
            rgtag = rgtag->next;
        }
        rgArray[n] = entry;
    }

    // Remove the existing RG lines
    sh = sam_hdr_del(sh, "RG", NULL, NULL);

    // add the new RG lines
    for (n=0; n<nrg; n++) {
        char *entry = strdup(rgArray[n]);
        addNewRG(sh, entry, "0", NULL, NULL, NULL);
        free(entry);

        // for each tag in barcodeArray
        for (i=1; i < barcodeArray->end; i++) {
            bc_details_t *bcd = barcodeArray->entries[i];

            char *entry = strdup(rgArray[n]);
            addNewRG(sh, entry, bcd->name, bcd->lib, bcd->sample, bcd->desc);
            free(entry);
        }
    }

    for (n=0; n<nrg; n++) {
        free(rgArray[n]);
    }
    free(rgArray);

    free(h->text);
    sam_hdr_rebuild(sh);
    h->text = strdup(sam_hdr_str(sh));
    h->l_text = sam_hdr_length(sh);
    sam_hdr_free(sh);
}

/*
 * Process one template
 */
static int processTemplate(va_t *template, va_t *barcodeArray, HashTable *barcodeHash, HashTable *tagHopHash, decode_opts_t *opts)
{
    char *name = NULL;
    char *bc_tag = NULL;
    char *qt_tag = NULL;
    char *newtag = NULL;
    short error_code = 0;
    size_t bc_len = 0;
    char stack_bc_tag[STACK_BC_LEN];
    char stack_qt_tag[STACK_BC_LEN];
    char stack_newtag[STACK_BC_LEN];

    // look for barcode tag
    for (int n=0; n < template->end; n++) {
        bam1_t *rec = template->entries[n];
        uint8_t *p = bam_aux_get(rec,opts->barcode_tag_name);
        if (p) {
            if (bc_tag) { // have we already found a tag?
                if (strcmp(bc_tag,bam_aux2Z(p)) != 0) {
                    fprintf(stderr,"Record %s has two different barcode tags: %s and %s\n",
                                   bam_get_qname(rec), bc_tag, bam_aux2Z(p));
                    if (bc_tag != stack_bc_tag) free(bc_tag);
                    if (qt_tag != stack_qt_tag) free(qt_tag);
                    return -1;
                }
            } else {
                char *bc = bam_aux2Z(p);
                bc_len = strlen(bc);
                if (bc_len < sizeof(stack_bc_tag)) {
                    bc_tag = stack_bc_tag;
                    memcpy(bc_tag, bc, bc_len + 1);
                } else {
                    bc_tag = strdup(bc);
                    if (!bc_tag) die("Out of memory");
                }
                p = bam_aux_get(rec,opts->quality_tag_name);
                if (p) {
                    char *qt = bam_aux2Z(p);
                    size_t qt_len = strlen(qt);
                    if (qt_len < sizeof(stack_qt_tag)) {
                        qt_tag = stack_qt_tag;
                        memcpy(qt_tag, qt, qt_len + 1);
                    } else {
                        qt_tag = strdup(bam_aux2Z(p));
                        if (!qt_tag) die("Out of memory");
                    }
                }
            }
        }
    }

    // if the convert_low_quality flag is set, then (potentially) change the tag
    if (bc_tag) {
        if (opts->convert_low_quality && qt_tag) {
            if (bc_len >= sizeof(stack_newtag)) {
                newtag = malloc(bc_len + 1);
                if (!newtag) die("Out of memory");
            } else {
                newtag = stack_newtag;
            }
            checkBarcodeQuality(newtag, bc_tag, bc_len, qt_tag, strlen(qt_tag), opts);
        } else {
            newtag = bc_tag;
        }
        // truncate to barcode lengths if necessary
        char stack_idx1[STACK_BC_LEN];
        char stack_idx2[STACK_BC_LEN];
        char *idx1 = stack_idx1, *idx2 = stack_idx2;
        split_index(bc_tag, bc_len, opts->dual_tag, &idx1, &idx2, sizeof(stack_idx1), sizeof(stack_idx2));
        if ( (strlen(idx1) > opts->idx1_len) || (strlen(idx2) > opts->idx2_len) ) {
            if (strlen(idx1) > opts->idx1_len) idx1[opts->idx1_len]=0;
            if (strlen(idx2) > opts->idx2_len) idx2[opts->idx2_len]=0;
            strcpy(newtag,idx1); 
            if (opts->idx2_len) strcat(newtag,INDEX_SEPARATOR);
            strcat(newtag,idx2);
        }
        if (idx1 != stack_idx1) free(idx1);
        if (idx2 != stack_idx2) free(idx2);
    }

    for (int n=0; n < template->end; n++) {
        bam1_t *rec = template->entries[n];
        if (newtag) {
            char stack_newrg[256];
            char *newrg = stack_newrg;
            if (n==0) name = findBarcodeName(newtag,barcodeArray, barcodeHash, tagHopHash, opts,!(rec->core.flag & BAM_FQCFAIL), n==0);
            makeNewTag(rec,"RG",name, &newrg, sizeof(stack_newrg));
            bam_aux_update_str(rec,"RG",strlen(newrg)+1, newrg);
            if (newrg != stack_newrg) free(newrg);
            if (opts->change_read_name) add_suffix(rec, name);
        }
    }

    if (newtag != bc_tag && newtag != stack_newtag) free(newtag);
    if (qt_tag && qt_tag != stack_qt_tag) free(qt_tag);
    if (bc_tag && bc_tag != stack_bc_tag) free(bc_tag);
    return error_code;
}

/*
 * Read records from a given iterator until the qname changes
 */
static va_t *loadTemplate(BAMit_t *bit, char *qname)
{
    va_t *recordSet = va_init(5,freeRecord);

    while (BAMit_hasnext(bit) && strcmp(bam_get_qname(BAMit_peek(bit)),qname) == 0) {
        bam1_t *rec = bam_init1();
        bam_copy1(rec,BAMit_next(bit));
        va_push(recordSet,rec);
    }

    return recordSet;
}

static int processTemplatesNoThreads(BAMit_t *bam_in, BAMit_t *bam_out, va_t *barcodeArray, HashTable *barcodeHash, HashTable *tagHopHash, decode_opts_t* opts)
{
    char qname[257] = {0};

    while (BAMit_hasnext(bam_in)) {
        bam1_t *rec = BAMit_peek(bam_in);
        va_t *template;
        memcpy(qname, bam_get_qname(rec), rec->core.l_qname);
        template = loadTemplate(bam_in, qname);
        if (processTemplate(template, barcodeArray, barcodeHash, tagHopHash, opts)) break;
        for (int n = 0; n < template->end; n++) {
            bam1_t *rec_n = template->entries[n];
            int r = sam_write1(bam_out->f, bam_out->h, rec_n);
            if (r < 0) {
                fprintf(stderr, "Could not write sequence\n");
                return -1;
            }
        }
        va_free(template);
    }
    return 0;
}

static void *decode_job(void *arg)
{
    decode_thread_data_t *job_data = (decode_thread_data_t *) arg;
    va_t template = { 0, 0, NULL, NULL };
    size_t start_rec = 0;

    job_data->result = -1;
    for (int i = 0; i < job_data->template_counts->end; i++) {
        template.end = template.max = job_data->template_counts->entries[i];
        template.entries = &job_data->record_set->entries[start_rec];
        start_rec += template.end;
        if (processTemplate(&template, job_data->barcode_array, job_data->barcodeHash, job_data->tagHopHash, job_data->opts)) goto fail;
    }
    assert(start_rec == job_data->nrec);

    job_data->result = 0;
 fail:
    return job_data;
}

static void output_job_results(BAMit_t *bam_out, decode_thread_data_t *job_data)
{
    if (job_data->result != 0) {
        die("Processing job failed to return a result\n");
    }

    // Write out result records
    for (int i = 0; i < job_data->nrec; i++) {
        bam1_t *rec = job_data->record_set->entries[i];
        int r = sam_write1(bam_out->f, bam_out->h, rec);
        if (r < 0) {
            die("Could not write sequence\n");
        }
    }
}

void accumulate_job_metrics(va_t *job_barcodes, HashTable *job_tag_hops, va_t *barcodeArray, HashTable *tagHopHash)
{
    // Accumulate metrics
    HashIter *iter = HashTableIterCreate();
    if (!iter) die("Out of memory");
    HashItem *hi;

    for (int i = 0; i < barcodeArray->end; i++) {
        bc_details_t *bc = barcodeArray->entries[i];
        bc_details_t *job_bc = job_barcodes->entries[i];
        bc->reads           += job_bc->reads;
        bc->pf_reads        += job_bc->pf_reads;
        bc->perfect         += job_bc->perfect;
        bc->pf_perfect      += job_bc->pf_perfect;
        bc->one_mismatch    += job_bc->one_mismatch;
        bc->pf_one_mismatch += job_bc->pf_one_mismatch;
    }

    // Accumulate tag hops
    while ((hi = HashTableIterNext(job_tag_hops, iter)) != NULL) {
        int added = -1;
        HashItem *hi2 = HashTableAdd(tagHopHash, hi->key, hi->key_len, hi->data, &added);
        if (!hi2) die("Out of memory");
        if (!added) {  // Already in there so accumulate results
            bc_details_t *job_bc = hi->data.p;
            bc_details_t *acc_bc = hi2->data.p;
            acc_bc->reads           += job_bc->reads;
            acc_bc->pf_reads        += job_bc->pf_reads;
            acc_bc->perfect         += job_bc->perfect;
            acc_bc->pf_perfect      += job_bc->pf_perfect;
            acc_bc->one_mismatch    += job_bc->one_mismatch;
            acc_bc->pf_one_mismatch += job_bc->pf_one_mismatch;
            free_taghop_bcd(job_bc);
        }
    }

    HashTableIterDestroy(iter);
}

va_t *copy_barcode_array(va_t *barcode_array)
{
    va_t *copy = va_init(barcode_array->end, NULL);
    bc_details_t *barcodes;

    if (barcode_array->end == 0) return copy;

    barcodes = calloc(barcode_array->end, sizeof(*barcodes));
    if (!barcodes) die("Out of memory");

    for (int i = 0; i < barcode_array->end; i++) {
        bc_details_t *bc = barcode_array->entries[i];
        barcodes[i].seq    = bc->seq;
        barcodes[i].idx1   = bc->idx1;
        barcodes[i].idx2   = bc->idx2;
        barcodes[i].name   = bc->name;
        barcodes[i].lib    = bc->lib;
        barcodes[i].sample = bc->sample;
        barcodes[i].desc   = bc->desc;
        va_push(copy, &barcodes[i]);
    }
    return copy;
}

void delete_barcode_array_copy(va_t *barcode_array)
{
    if (barcode_array->end > 0) free(barcode_array->entries[0]);
    va_free(barcode_array);
}

static void job_free(decode_thread_data_t *job_data)
{
    va_free(job_data->record_set);
    ia_free(job_data->template_counts);
    delete_barcode_array_copy(job_data->barcode_array);
    HashTableDestroy(job_data->tagHopHash, 0);
    free(job_data);
}

static decode_thread_data_t *init_job(va_t *barcode_array, HashTable *barcodeHash, decode_opts_t *opts)
{
    decode_thread_data_t *job_data = calloc(1, sizeof(*job_data));
    if (!job_data) die("Out of memory\n");

    job_data->barcode_array = copy_barcode_array(barcode_array);
    job_data->barcodes = job_data->barcode_array->end > 0 ? job_data->barcode_array->entries[0] : NULL;
    job_data->tagHopHash = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS);
    if (!job_data->tagHopHash) die("Out of memory");
    job_data->barcodeHash = barcodeHash;
    job_data->opts = opts;
    job_data->result = -1;
    job_data->next = NULL;

    return job_data;
}

static int processTemplatesThreads(hts_tpool *pool, BAMit_t *bam_in, BAMit_t *bam_out, va_t *barcodeArray, HashTable *barcodeHash, HashTable *tagHopHash, decode_opts_t* opts)
{
    hts_tpool_result *job_result = NULL;
    hts_tpool_process *queue = hts_tpool_process_init(pool, 2 * opts->nthreads, 0);
    decode_thread_data_t *job_freelist = NULL;
    decode_thread_data_t *job_data = init_job(barcodeArray, barcodeHash, opts);
    char qname[257] = { 0 };

    if (!queue) {
        perror("hts_tpool_process_init");
        return -1;
    }

    job_data->record_set = va_init(TEMPLATES_PER_JOB * 2, freeRecord);
    job_data->template_counts = ia_init(TEMPLATES_PER_JOB);
    job_data->nrec = 0;

    while (BAMit_hasnext(bam_in)) {
        bam1_t *rec = BAMit_peek(bam_in);
        int rec_count = 0;
        memcpy(qname, bam_get_qname(rec), rec->core.l_qname);
        while (BAMit_hasnext(bam_in) && strcmp(bam_get_qname(BAMit_peek(bam_in)),qname) == 0) {
            if (job_data->nrec < job_data->record_set->end) {
                bam_copy1(job_data->record_set->entries[job_data->nrec], BAMit_next(bam_in));
            } else {
                bam1_t *rec = bam_init1();
                if (!rec) { die("Out of memory"); }
                bam_copy1(rec, BAMit_next(bam_in));
                va_push(job_data->record_set, rec);
            }
            job_data->nrec++;
            rec_count++;
        }
        ia_push(job_data->template_counts, rec_count);

        if (job_data->template_counts->end == TEMPLATES_PER_JOB) {
            while (job_data != NULL) {
                int blk = hts_tpool_dispatch2(pool, queue, decode_job, job_data, 1);
                if (!blk) {
                    job_data = NULL;
                } else if (errno != EAGAIN) {
                    die("Thread pool dispatch failed");
                }

                if (blk) {
                    job_result = hts_tpool_next_result_wait(queue);
                    if (!job_result) {
                        die("Failed to get processing job result");
                    }
                } else {
                    job_result = hts_tpool_next_result(queue);
                }

                if (job_result != NULL) {
                    decode_thread_data_t *finished_job = hts_tpool_result_data(job_result);
                    output_job_results(bam_out, finished_job);
                    finished_job->next = job_freelist;
                    job_freelist = finished_job;
                    hts_tpool_delete_result(job_result, 0);
                }
            }

            if (job_freelist != NULL) {
                job_data = job_freelist;
                job_freelist = job_data->next;
                job_data->template_counts->end = 0;
            }  else {
                job_data = init_job(barcodeArray, barcodeHash, opts);
                job_data->record_set = va_init(TEMPLATES_PER_JOB * 2, freeRecord);
                job_data->template_counts = ia_init(TEMPLATES_PER_JOB);
            }
            job_data->nrec = 0;
        }
    }

    if (job_data->template_counts->end > 0) {
        // Deal with left-over items
        if (hts_tpool_dispatch(pool, queue, decode_job, job_data) < 0) {
            die("Thread pool dispatch failed");
        }
    } else {
        // Not used, put back on free list
        job_data->next = job_freelist;
        job_freelist = job_data;
    }

    while (!hts_tpool_process_empty(queue)) {
        job_result = hts_tpool_next_result_wait(queue);
        if (!job_result) {
            die("Failed to get processing job result");
        }
        decode_thread_data_t *finished_job = hts_tpool_result_data(job_result);
        output_job_results(bam_out, finished_job);
        finished_job->next = job_freelist;
        job_freelist = finished_job;
        hts_tpool_delete_result(job_result, 0);
    }

    while (job_freelist != NULL) {
        decode_thread_data_t *next = job_freelist->next;
        accumulate_job_metrics(job_freelist->barcode_array, job_freelist->tagHopHash, barcodeArray, tagHopHash);
        job_free(job_freelist);
        job_freelist = next;
    }

    hts_tpool_process_destroy(queue);

    return 0;
}

/*
 * Make a hash table from barcodeArray
 */
HashTable *make_barcode_hash(va_t *barcodeArray)
{
    HashTable *barcodeHash = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS);
    if (!barcodeHash) die("out of memory");
    for (int n=0; n < barcodeArray->end; n++) {
        bc_details_t *bcd = barcodeArray->entries[n];
        HashData hd;
        hd.i = n;
        if (HashTableAdd(barcodeHash, bcd->seq, 0, hd, NULL) == NULL) {
            die("Out of memory");
        }
    }
    return barcodeHash;
}

/*
 * Returns the length of the longest name
 */
size_t find_longest_barcode_name(va_t *barcodeArray)
{
    size_t longest = 0;
    for (int n=0; n < barcodeArray->end; n++) {
        bc_details_t *bcd = barcodeArray->entries[n];
        size_t len = bcd->name ? strlen(bcd->name) : 0;
        if (len > longest) longest = len;
    }
    return longest;
}

int get_barcode_metadata(va_t *barcodeArray, int idx,
                         const char **name_out, const char **lib_out,
                         const char **sample_out, const char **desc_out) {
    bc_details_t *barcode;
    if (!barcodeArray || idx >= barcodeArray->end) return -1;
    barcode = barcodeArray->entries[idx];
    if (name_out)   *name_out   = barcode->name;
    if (lib_out)    *lib_out    = barcode->lib;
    if (sample_out) *sample_out = barcode->sample;
    if (desc_out)   *desc_out   = barcode->desc;
    return 0;
}

/*
 * Main code
 */
static int decode(decode_opts_t* opts)
{
    int retcode = 1;
    BAMit_t *bam_in = NULL;
    BAMit_t *bam_out = NULL;
    va_t *barcodeArray = NULL;
    HashTable *tagHopHash = NULL;
    HashTable *barcodeHash = NULL;
    htsThreadPool hts_threads = { NULL, 0 };

    while (1) {
        if (opts->nthreads > 1) {
            hts_threads.pool = hts_tpool_init(opts->nthreads);
            if (!hts_threads.pool) {
                fprintf(stderr, "Couldn't set up thread pool\n");
                break;
            }
        }

        /*
         * Read the barcode (tags) file 
         */
        barcodeArray = loadBarcodeFile(opts);
        if (!barcodeArray) break;

        // create hash from barcodeArray
        barcodeHash = make_barcode_hash(barcodeArray);

        tagHopHash = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS);

        /*
         * Open input fnd output BAM files
         */
        bam_in = BAMit_open(opts->input_name, 'r', opts->input_fmt, 0, hts_threads.pool ? &hts_threads : NULL);
        if (!bam_in) break;
        bam_out = BAMit_open(opts->output_name, 'w', opts->output_fmt, opts->compression_level, hts_threads.pool ? &hts_threads : NULL);
        if (!bam_out) break;
        // copy input to output header
        bam_hdr_destroy(bam_out->h); bam_out->h = bam_hdr_dup(bam_in->h);

        // Change header by adding PG and RG lines
        changeHeader(barcodeArray, bam_out->h, opts->argv_list);
        if (sam_hdr_write(bam_out->f, bam_out->h) != 0) {
            fprintf(stderr, "Could not write output file header\n");
            break;
        }

        // Read and process each template in the input BAM
        if (opts->nthreads < 2) {
            if (processTemplatesNoThreads(bam_in, bam_out, barcodeArray, barcodeHash, tagHopHash, opts) < 0) break;
        } else {
            if (processTemplatesThreads(hts_threads.pool, bam_in, bam_out, barcodeArray, barcodeHash, tagHopHash, opts) < 0) break;
        }

        if (BAMit_hasnext(bam_in)) break;   // we must has exited the above loop early

        /*
         * And finally.....the metrics
         */
        if (opts->metrics_name) {
            if (writeMetrics(barcodeArray, tagHopHash, opts) != 0) break;
        }
                
        retcode = 0;
        break;
    }

    // tidy up after us
    HashIter *iter = HashTableIterCreate();
    if (iter) {
        HashItem *hi;
        while ((hi = HashTableIterNext(tagHopHash, iter)) != NULL) {
            free_taghop_bcd(hi->data.p);
        }
        HashTableIterDestroy(iter);
    }
    va_free(barcodeArray);
    HashTableDestroy(barcodeHash, 0);
    HashTableDestroy(tagHopHash, 0);
    BAMit_free(bam_in);
    BAMit_free(bam_out);
    if (hts_threads.pool) hts_tpool_destroy(hts_threads.pool);

    return retcode;
}

/*
 * called from bambi to perform index decoding
 *
 * Parse the command line arguments, then call the main decode function
 *
 * returns 0 on success, 1 if there was a problem
 */
int main_decode(int argc, char *argv[])
{
    int ret = 1;

    decode_opts_t* opts = parse_args(argc, argv);
    if (opts) {
        ret = decode(opts);
    }
    decode_free_opts(opts);
    return ret;
}
