/*  adapter.c -- index adapter subcommand.

    Copyright (C) 2019 Genome Research Ltd.

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
#include <htslib/thread_pool.h>
#include <htslib/hfile.h>
#include <cram/sam_header.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>

#include "bamit.h"
#include "parse_bam.h"
#include "adapters.h"
#include "hash_table.h"

#define DEBUG 1

#define xstr(s) str(s)
#define str(s) #s

#define min(a, b) ( (a<=b) ? a : b )

#define DEFAULT_MIN_SCORE 16
#define DEFAULT_MIN_FRAC .75
#define DEFAULT_MIN_PFRAC .8

#define TEMPLATES_PER_JOB 5000

#define SEEDLEN 12
#define MAXSEEDDIFF 2

uint8_t S[256];
uint64_t SEEDMASK;

// Metrics are stored in a global hash table keyed on RG tag
typedef struct {
    int total_fwd;
    int total_rev;
    int contam_fwd;
    int contam_rev;
    int hist_len;
    int *hist_fwd;
    int *hist_rev;
    HashTable *adapter; // key is adapter name, value is adapter_metrics_t
} rg_metrics_t;

typedef struct {
    int fwd, rev;
} adapter_metrics_t;

static HashTable *metrics = NULL;
static pthread_mutex_t metrics_lock = PTHREAD_MUTEX_INITIALIZER;

/*
 * structure to hold options
 */
typedef struct {
    char *input_name;
    char *output_name;
    char *metrics_name;
    va_t *adapterArray;
    bool verbose;
    char *argv_list;
    char *input_fmt;
    char *output_fmt;
    char compression_level;
    int nthreads;
    int minscore;
    double minfrac;
    double minpfrac;
} adapter_opts_t;

static void freeAdapter(void *r)
{
    adapter_t *a = (adapter_t *)r;
    if (!a) return;
    free(a->name);
    free(a->fwd);
    free(a->rev);
    free(a);
}

/*
 * trim whitespace from string
 */
static void trim(char *str)
{
    char *end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;
    end[1] = 0;
}

/*
 * Read the adapter array from a FASTA file
 */
static va_t *loadAdapterFile(char *fname)
{
    char buf[1024];
    size_t n = 1024;
    va_t *a = va_init(500,freeAdapter);

    hFILE *f = hopen(fname, "r");
    if (!f) die("Can't open adapter file '%s'", fname);

    while (hgets(buf,n,f) > 0) {
        adapter_t *adapter = smalloc(sizeof(adapter_t));
        if (*buf != '>') die("problem reading adapter file '%s'\nThis does not look like a FASTA file", fname);
        trim(buf); adapter->name = strdup(buf+1);
        if (hgets(buf,n,f) <= 0) die("Incomplete adapter file");
        trim(buf); adapter->fwd = strdup(buf);
        adapter->rev = strdup(buf);
        rev_comp_seq(adapter->rev);
        adapter->offset = 0;
        va_push(a,adapter);
    }

    if (hclose(f)) die("Can't close adapter file '%s'", fname);
    return a;
}

/*
 * Create Kmer fragments for each adapter
 */
static void fragmentAdapters(va_t *adapters, int minscore)
{
    int N = adapters->end;
    for (int n=0; n < N; n++) {
        adapter_t *a = adapters->entries[n];
        for (int pos = 1; pos < strlen(a->fwd) - minscore; pos++) {
            adapter_t *newAdapter = smalloc(sizeof(adapter_t));
            newAdapter->name = strdup(a->name);
            newAdapter->fwd = strdup(a->fwd+pos);
            newAdapter->rev = strdup(a->rev+pos);
            newAdapter->offset = pos;
            va_push(adapters, newAdapter);
        }
    }
}

/*
 * Initialise various seed values
 */
static void initSeedValues(void)
{
    memset(S, 4, sizeof(S));
    S['a'] = S['A'] = 0;
    S['c'] = S['C'] = 1;
    S['g'] = S['G'] = 2;
    S['t'] = S['T'] = 3;

    SEEDMASK =
        (1ull << 0) |
        (1ull << 3) |
        (1ull << 6) |
        (1ull << 9) |
        (1ull << 12) |
        (1ull << 15) |
        (1ull << 18) |
        (1ull << 21) |
        (1ull << 24) |
        (1ull << 27) |
        (1ull << 30) |
        (1ull << 33) |
        (1ull << 36) |
        (1ull << 39) |
        (1ull << 42) |
        (1ull << 45) |
        (1ull << 48) |
        (1ull << 51) |
        (1ull << 54) |
        (1ull << 57) |
        (1ull << 60) |
        (1ull << 63);

}

static inline int seedDiff(uint64_t s1, uint64_t s2)
{
    uint64_t dif = (s1 ^ s2);
    dif = (dif | (dif >> 1) | (dif >> 2)) & SEEDMASK;
    return __builtin_popcountll(dif);
}

/*
 * Calculate seed for a string
 */
static uint64_t calcSeed(char *seq)
{
    uint64_t seed = 0;
    if (strlen(seq) < SEEDLEN) return seed;
    for (int n=0; n < SEEDLEN; n++) {
        seed = seed << 3;
        seed |= S[(int)*(seq++)];
    }
    return seed;
}


/*
 * Calculate seed for each adapter fragment
 */
static void calcAdapterSeed(va_t *adapters)
{
    for (int n=0; n < adapters->end; n++) {
        adapter_t *a = adapters->entries[n];
        a->fwd_seed = calcSeed(a->fwd);
        a->rev_seed = calcSeed(a->rev);
    }
}


/*
 * compare two adapters a,b
 * return -1 if a < b
 *         0 if a == b
 *        +1 if a > b
 */
static int compareAdapters(adapter_t *a, adapter_t *b)
{
    if (!b) return 1;
    if (a->score == b->score) {
        if (a->offset == b->offset) return 0;
        if (a->offset < b->offset) return 1;
        return -1;
    }
    if (a->score > b->score) return 1;
    return -1;
}

static adapter_opts_t *adapter_init_opts(int argc, char **argv)
{
    adapter_opts_t *opts = calloc(1, sizeof(*opts));
    if (!opts) die("Out of memory");
    opts->argv_list = stringify_argv(argc, argv);
    if (opts->argv_list[strlen(opts->argv_list)-1] == ' ') opts->argv_list[strlen(opts->argv_list)-1] = 0;
    opts->adapterArray = NULL;

    // set defaults
    opts->minscore = DEFAULT_MIN_SCORE;
    opts->minfrac = DEFAULT_MIN_FRAC;
    opts->minpfrac = DEFAULT_MIN_PFRAC;
    opts->verbose = false;
    return opts;
}

static void adapter_free_opts(adapter_opts_t *opts)
{
    if (!opts) return;
    free(opts->input_name);
    free(opts->output_name);
    free(opts->argv_list);
    free(opts->input_fmt);
    free(opts->output_fmt);
    free(opts->metrics_name);
    va_free(opts->adapterArray);
    free(opts);
}

// Data for thread pool jobs
typedef struct adapter_thread_data_t {
    va_t *record_set;                   // records to process
    ia_t *template_counts;              // records in each template
    va_t *adapter_array;                // job-local copy of adapters array
    adapter_t *adapter;             // memory for job-local adapter
    adapter_opts_t *opts;                       // pointer to shared opts
    int nrec;                           // number of live records in record_set
    int result;                         // job result, 0 = success
    struct adapter_thread_data_t *next;  // for free list
    BAMit_t *bam_out;
} adapter_thread_data_t;

static void freeRecord(void *r) { bam_destroy1((bam1_t *)r); }

//
// resize histogram arrays
//
static void resize_hist(rg_metrics_t *rgm, int newsize)
{
    newsize++;
    rgm->hist_fwd = realloc(rgm->hist_fwd, sizeof(*rgm->hist_fwd)*newsize);
    if (!rgm->hist_fwd) die("Out of memory in resize_hist(%d)", newsize);

    rgm->hist_rev = realloc(rgm->hist_rev, sizeof(*rgm->hist_rev)*newsize);
    if (!rgm->hist_rev) die("Out of memory in resize_hist(%d)", newsize);

    // initialise new values
    for (int n = rgm->hist_len; n < newsize; n++) {
        rgm->hist_fwd[n] = 0;
        rgm->hist_rev[n] = 0;
    }
    rgm->hist_len = newsize;
}

static void updateMetrics(HashTable *m, bam1_t *rec, adapter_t *adapter)
{
    HashItem *hi;
    rg_metrics_t *rgm = NULL;
    uint8_t *s = bam_aux_get(rec, "RG");
    if (!s) return;

    char *rg = bam_aux2Z(s);
    hi = HashTableSearch(m, rg, strlen(rg));
    if (!hi) {
        HashData hd;
        rgm = calloc(1, sizeof(*rgm));
        hd.p = rgm;
        HashTableAdd(m, rg, strlen(rg), hd, NULL);
    } else {
        rgm = (rg_metrics_t *)hi->data.p;
    }

    if (bam_is_rev(rec)) {
        rgm->total_rev++;
        if (adapter) {
            rgm->contam_rev++;
            if (adapter->seqstart >= rgm->hist_len) { resize_hist(rgm, adapter->seqstart); }
            rgm->hist_rev[adapter->seqstart]++;
        }
    } else {
        rgm->total_fwd++;
        if (adapter) {
            rgm->contam_fwd++;
            if (adapter->seqstart >= rgm->hist_len) { resize_hist(rgm, adapter->seqstart); }
            rgm->hist_fwd[adapter->seqstart]++;
        }
    }
}

/*
 * display usage information
 */
static void usage(FILE *write_to)
{
    fprintf(write_to,
"Usage: bambi adapter [options] filename\n"
"\n"
"Options:\n"
"  -o   --output                        output file [default: stdout]\n"
"  -v   --verbose                       verbose output\n"
"       --metrics-file                  metrics written to this file\n"
"  -a   --adapter-file                  use file of adapters instead of built-in list\n"
"                                       The file must be in FASTA format\n"
"       --input-fmt                     format of input file [sam/bam/cram]\n"
"       --output-fmt                    format of output file [sam/bam/cram]\n"
"       --compression-level             Compression level of output file [0..9]\n"
"  -t   --threads                       number of threads to use [default: 1]\n"
);
}

/*
 * Takes the command line options and turns them into something we can understand
 */
static adapter_opts_t* parse_args(int argc, char *argv[])
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char* optstring = "a:i:o:vb:t:";

    static const struct option lopts[] = {
        { "input",                      1, 0, 'i' },
        { "output",                     1, 0, 'o' },
        { "verbose",                    0, 0, 'v' },
        { "metrics-file",               1, 0, 0 },
        { "adapter-file",               1, 0, 'a' },
        { "input-fmt",                  1, 0, 0 },
        { "output-fmt",                 1, 0, 0 },
        { "compression-level",          1, 0, 0 },
        { "threads",                    1, 0, 't' },
        { NULL, 0, NULL, 0 }
    };

    adapter_opts_t* opts = adapter_init_opts(argc+1, argv-1);
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
        case 't':   opts->nthreads = atoi(optarg);
                    break;
        case 'a':   opts->adapterArray = loadAdapterFile(optarg);
                    break;
        case 0:     arg = lopts[option_index].name;
                         if (strcmp(arg, "metrics-file") == 0)               opts->metrics_name = strdup(optarg);
                    else if (strcmp(arg, "input-fmt") == 0)                  opts->input_fmt = strdup(optarg);
                    else if (strcmp(arg, "output-fmt") == 0)                 opts->output_fmt = strdup(optarg);
                    else if (strcmp(arg, "compression-level") == 0)          opts->compression_level = *optarg;
                    else {
                        printf("\nUnknown option: %s\n\n", arg); 
                        usage(stdout); adapter_free_opts(opts);
                        return NULL;
                    }
                    break;
        default:    printf("Unknown option: '%c'\n", opt);
            /* else fall-through */
        case '?':   usage(stdout); adapter_free_opts(opts); return NULL;
        }
    }

    argc -= optind;
    argv += optind;

    if (argc > 0) opts->input_name = strdup(argv[0]);
    optind = 0;

    // input defaults to stdin
    if (!opts->input_name) opts->input_name = strdup("-");

    // output defaults to stdout
    if (!opts->output_name) opts->output_name = strdup("-");

    // use built-in adapters if adapter file not specified
    if (!opts->adapterArray) {
        opts->adapterArray = va_init(500,freeAdapter);
        adapter_t *src = DefaultAdapterArray;
        while (src->name) {
            adapter_t *dest = smalloc(sizeof(adapter_t));
            dest->name = strdup(src->name);
            dest->fwd = strdup(src->fwd);
            dest->rev = strdup(src->rev);
            dest->offset = 0;
            va_push(opts->adapterArray, dest);
            src++;
        }
    }

    fragmentAdapters(opts->adapterArray, opts->minscore);
    initSeedValues();
    calcAdapterSeed(opts->adapterArray);

    return opts;
}

static void initAdapterScore(adapter_t *a)
{
    a->len = strlen(a->fwd);
    a->score = 0; a->begin = 0; a->end = a->len - 1;
    a->frac = 0.0; a->pfrac = 0.0;
}

static double kmerPoisson(bam1_t *rec, adapter_t *adapter)
{
    uint64_t fA=0, fC=0, fG=0, fT=0;

    char *aseq = adapter->revmatch ? adapter->rev : adapter->fwd;
    int matchlen = rec->core.l_qseq - adapter->seqstart;
    int adplen = adapter->len - adapter->begin;
    int comlen = min(matchlen,adplen);

    for (int n=0; n < comlen; n++) {
        char b = aseq[n + adapter->begin];
        switch (b) {
            case 'a': case 'A': fA++; break;
            case 'c': case 'C': fC++; break;
            case 'g': case 'G': fG++; break;
            case 't': case 'T': fT++; break;
            default: switch (random() % 4) {
                        case 0: fA++; break;
                        case 1: fC++; break;
                        case 2: fG++; break;
                        case 3: fT++; break;
                     }
        }
    }

    uint64_t L = 3000000000ull;
    uint64_t K = fA + fC + fG + fT;
    double lambda = (L+1-K) * pow(0.25,fA) * pow(0.25,fC) * pow(0.25,fG) * pow(0.25,fT);
    return 1 / exp(lambda);
}

static void dumpAdapterResult(bam1_t *rec, adapter_t *a)
{
    char *seq = (char *)get_read(rec);
    char *aseq = (a->revmatch ? a->rev : a->fwd);
    int len = min(strlen(seq), a->len);

    fprintf(stderr,"sco: %d  off: %d  beg: %d  end: %d  len:%d  mat:%d  ", a->score, a->offset, a->begin, a->end, a->end-a->begin, a->begin+a->seqstart);
    fprintf(stderr,"af:%f  frac:%f\n", a->pfrac, a->frac);
    fprintf(stderr,"%s\n", bam_get_qname(rec));
    fprintf(stderr,"%s\n", a->name);
    fprintf(stderr,"%.*s\n", len, seq + a->seqstart);
    fprintf(stderr,"%s\n", aseq);
    for (int n=0; n < len; n++) {
        fprintf(stderr, "%c", (seq[n+a->seqstart]==aseq[n]) ? '+' : '-');
    }
    fprintf(stderr, "\n");
    for (int n=0; n < len; n++) {
        fprintf(stderr, "%c", (n<a->begin || n>=a->end) ? ' ' : '*');
    }
    fprintf(stderr, "\n\n");
    free(seq);
}

#define SCORE_MATCH 1
#define PEN_MISMATCH -2

static void calcAdapterScore(char *seq, int seqstart, char *aseq, adapter_t *a, adapter_opts_t *opts, bool rev)
{
    uint64_t w = 0;
    for (int n = 0; n < SEEDLEN; n++) {
        w = w << 3;
        w |= S[(int)seq[seqstart+n]];
    }

    // if seeds do not match, get out quick
    if (seedDiff(w, (rev ? a->rev_seed : a->fwd_seed)) > MAXSEEDDIFF) return;

    int comlen = min(strlen(seq+seqstart), a->len);
    int score = 0, maxscore = 0;
    int currstart = 0, maxstart = 0, maxend = 0;

    for (int i = 0; i < comlen; i++) {
        if (seq[seqstart+i] == aseq[i]) score += SCORE_MATCH;
        else                            score += PEN_MISMATCH;

        if (score < 0) {
            score = 0;
            currstart = i + 1;
        } else {
            if (score > maxscore) {
                maxscore = score;
                maxstart = currstart;
                maxend = i + 1;
            }
        }
    }

    if (maxscore > a->score) {
        a->score = maxscore;
        a->begin = maxstart;
        a->end = maxend;
        a->revmatch = rev;
        a->seqstart = seqstart;
        a->frac = (double) (a->end - a->begin) / (double) (a->len + a->offset);
        a->pfrac = (double) (a->end - a->begin) / (double) (comlen + a->offset);
        // doesn't matter how good the score is, if it fails the other criteria
        if (a->frac < opts->minfrac || a->pfrac < opts->minpfrac) a->score = 0;
    }
}

static void matchAdapter(char *seq, adapter_t *adapter, adapter_opts_t *opts, bool isReverse)
{
    initAdapterScore(adapter);
    int seqStart = 0, seqEnd = strlen(seq) - 1 - opts->minscore;

    if (seqEnd < 0) seqEnd = 0;

    for (seqStart=0; seqStart < seqEnd; seqStart++) {
        calcAdapterScore(seq, seqStart, adapter->fwd, adapter, opts, false);
        calcAdapterScore(seq, seqStart, adapter->rev, adapter, opts, true);
    }
}

/*
 * find the adapter with the best match.
 */
static adapter_t *findBestMatch(bam1_t *rec, va_t *adapterArray, adapter_opts_t *opts)
{
    adapter_t *best_match = NULL;
    uint8_t *seq = get_read(rec);

    for (int n=0; n < adapterArray->end; n++) {
        adapter_t *adapter = adapterArray->entries[n];
        matchAdapter((char *)seq, adapter, opts, bam_is_rev(rec));
//dumpAdapterResult(rec,adapter);
        if (compareAdapters(adapter, best_match) > 0) {
            best_match = adapter;
        }
    }

    free(seq);
    return best_match;
}

static void updateRecord(bam1_t *rec, adapter_t *adapter)
{
    uint32_t clip = rec->core.l_qseq - adapter->seqstart + adapter->offset;
    float randconf = kmerPoisson(rec, adapter);
    bam_aux_append(rec, "aa", 'Z', strlen(adapter->name)+1, (uint8_t *) adapter->name);
    bam_aux_update_float(rec, "af", adapter->pfrac);
    bam_aux_update_float(rec, "ar", randconf);
    bam_aux_append(rec, "as", 'i', 4, (uint8_t *)&clip);
#if DEBUG
    uint32_t score = adapter->score;
    bam_aux_append(rec, "sc", 'i', 4, (uint8_t *)&score);
#endif
}

/*
 * Rewrite BAM header (with new PG line)
 */ 
static void changeHeader(bam_hdr_t *h, char *argv_list)
{
    SAM_hdr *sh = sam_hdr_parse_(h->text, h->l_text);

    sam_hdr_add_PG(sh, "bambi", "VN", bambi_version(), "CL", argv_list, NULL);

    free(h->text);
    sam_hdr_rebuild(sh);
    h->text = strdup(sam_hdr_str(sh));
    h->l_text = sam_hdr_length(sh);
    sam_hdr_free(sh);
}

/*
 * Check for overlap
 */
#define minoverlap 32
#define mismatchrate 0.1
#define almax 12

static void checkOverlap(va_t *template, bool verbose)
{
    bam1_t *rec1 = NULL;
    bam1_t *rec2 = NULL;
    char *seq1 = NULL;
    char *seq2 = NULL;

    for (int n=0; n < template->end; n++) {
        bam1_t *rec = template->entries[n];
        if (BAM_FREAD1 & rec->core.flag) rec1 = rec;
        if (BAM_FREAD2 & rec->core.flag) rec2 = rec;
    }

    // make sure we have a paired read
    if (rec1 == rec2) return;
    if (! (rec1 && rec2)) return;

    seq1 = (char *)get_read(rec1);
    seq2 = (char *)get_read(rec2);
    rev_comp_seq(seq2);

    uint64_t lseq1 = strlen(seq1);
    uint64_t lseq2 = strlen(seq2);
    uint64_t matchpos = lseq1;

    char *qe = seq1;
    char *q = qe + lseq1;

    do {
        uint64_t end1 = matchpos;
        uint64_t end2 = lseq2;
        uint64_t seedp1 = end1;
        uint64_t seedp2 = end2;
        uint64_t restoverlap = min(seedp1, seedp2);

        if (restoverlap >= minoverlap) {
            char *check1 = seq1 + seedp1 - restoverlap;
            char *check1e = check1 + restoverlap;
            char *check2 = seq2 + seedp2 - restoverlap;
            uint64_t maxmis = restoverlap * mismatchrate;
            uint64_t nummis = 0;

//fprintf(stderr,"%s\n", check1);
//fprintf(stderr,"%s\n", check2);
            while ( nummis <= maxmis && check1 != check1e )
                if (*check1++ != *check2++) nummis++;

//fprintf(stderr,"numis: %ld\n", nummis);
            if (check1 == check1e && nummis <= maxmis) {
                uint64_t al1 = (lseq1-end1);        // adapter length on read 1
                uint64_t al2 = end2 - restoverlap;  // adapter length on read 2
                uint64_t alcmp = min(almax, min(al1, al2));     // number of bases compared
                char *ap1 = seq1 + end1;        // adapter range on read 1
                char *ap1e = ap1 + alcmp;
                char *ap2 = seq2 + al2;
                unsigned int aldiff = 0;    // mismatches on adapter

                while (ap1 != ap1e) {
//fprintf(stderr,"%c %c\n", *ap1, *ap2);
                    if (*ap1++ != complement_base(*(--ap2))) aldiff++;
                }

                if (ap1 == ap1e && !aldiff) {
                    if (verbose) {
                        fprintf(stderr,"mismatchrate= %d / %d = %f", (int)nummis, (int)restoverlap, (double)nummis/restoverlap);
                        fprintf(stderr," al0=%d", (int)al1);
                        fprintf(stderr," al1=%d", (int)al2);
                        fprintf(stderr," aldiff=%d", aldiff);
                        fprintf(stderr," alcmp=%d", (int)alcmp);
                        fprintf(stderr,"\n");
                        fprintf(stderr,"[V2] assumed adapter on read 1 [%d] %.*s\n", (int)al1, (int)al1, seq1+end1);
                        fprintf(stderr,"[V2] assumed adapter on read 2 [%d] %.*s\n", (int)al2, (int)al2, seq2);
                    }
                    uint64_t ah = 1;
                    bam_aux_update_int(rec1, "ah", ah);
                    bam_aux_update_int(rec2, "ah", ah);
                    bam_aux_update_int(rec1, "a3", al1);
                    bam_aux_update_int(rec2, "a3", al2);
                }
            }
        }
        matchpos--;
        q--;
    } while (q != qe);

    free(seq2); free(seq1);
}

/*
 * Process one template
 */
static int processTemplate(va_t *template, va_t *adapterArray, adapter_opts_t *opts)
{
    int error_code = 0;
    adapter_t *adapter = NULL;

    for (int n=0; n < template->end; n++) {
        bam1_t *rec = template->entries[n];
        adapter = findBestMatch(rec, adapterArray, opts);
        if (adapter && adapter->score >= opts->minscore && adapter->frac >= opts->minfrac && adapter->pfrac >= opts->minpfrac) {
            updateRecord(rec, adapter);
            if (metrics) {
                pthread_mutex_lock(&metrics_lock);
                updateMetrics(metrics, rec, adapter);
                pthread_mutex_unlock(&metrics_lock);
            }
            if (opts->verbose) dumpAdapterResult(rec, adapter);
        } else {
            if (metrics) {
                pthread_mutex_lock(&metrics_lock);
                updateMetrics(metrics, rec, NULL);
                pthread_mutex_unlock(&metrics_lock);
            }
        }
    }

    checkOverlap(template, opts->verbose);

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

static int processTemplatesNoThreads(BAMit_t *bam_in, BAMit_t *bam_out, va_t *adapterArray, adapter_opts_t* opts)
{
    char qname[257] = {0};
    adapter_t *adapter = NULL;

    while (BAMit_hasnext(bam_in)) {
        bam1_t *rec = BAMit_peek(bam_in);
        va_t *template;
        memcpy(qname, bam_get_qname(rec), rec->core.l_qname);
        template = loadTemplate(bam_in, qname);

        // Look for adapters in the list
        for (int n=0; n < template->end; n++) {
            rec = template->entries[n];
            adapter = findBestMatch(rec, adapterArray, opts);
            if (adapter && adapter->score >= opts->minscore && adapter->frac >= opts->minfrac && adapter->pfrac >= opts->minpfrac) {
                updateRecord(rec, adapter);
                if (metrics) updateMetrics(metrics, rec, adapter);
                if (opts->verbose) dumpAdapterResult(rec, adapter);
            } else {
                if (metrics) updateMetrics(metrics, rec, NULL);
            }
        }
        // check for overlapping adapters
        checkOverlap(template, opts->verbose);

        // write records
        for (int n=0; n < template->end; n++) {
            rec = template->entries[n];
            int r = sam_write1(bam_out->f, bam_out->h, rec);
            if (r < 0) die("Could not write bam record\n");
        }
        va_free(template);
    }

    return 0;
}

static void *adapter_job(void *arg)
{
    adapter_thread_data_t *job_data = (adapter_thread_data_t *) arg;
    va_t template = { 0, 0, NULL, NULL };
    size_t start_rec = 0;

    job_data->result = -1;
    for (int i = 0; i < job_data->template_counts->end; i++) {
        template.end = template.max = job_data->template_counts->entries[i];
        template.entries = &job_data->record_set->entries[start_rec];
        start_rec += template.end;
        if (processTemplate(&template, job_data->adapter_array, job_data->opts)) goto fail;
    }
    assert(start_rec == job_data->nrec);

    job_data->result = 0;
 fail:
    return job_data;
}

static void output_job_results(BAMit_t *bam_out, adapter_thread_data_t *job_data)
{
    if (job_data->result != 0) {
        die("Processing job failed to return a result\n");
    }

    // Write out result records
    for (int i = 0; i < job_data->nrec; i++) {
        bam1_t *rec = job_data->record_set->entries[i];
        int r = sam_write1(bam_out->f, bam_out->h, rec);
        if (r < 0) die("Could not write sequence\n");
    }

}

static va_t *copy_adapter_array(va_t *adapter_array)
{
    va_t *dest = va_init(adapter_array->end, freeAdapter);

    for (int n = 0; n < adapter_array->end; n++) {
        adapter_t *src = adapter_array->entries[n];
        adapter_t *d = smalloc(sizeof(adapter_t));
        memset(d, 0, sizeof(adapter_t));
        d->name = strdup(src->name);
        d->fwd = strdup(src->fwd);
        d->rev = strdup(src->rev);
        d->offset = src->offset;
        d->fwd_seed = src->fwd_seed;
        d->rev_seed = src->rev_seed;
        va_push(dest,d);
    }
    return dest;
}

static void job_free(adapter_thread_data_t *job_data)
{
    va_free(job_data->record_set);
    ia_free(job_data->template_counts);
    va_free(job_data->adapter_array);
    free(job_data);
}

static adapter_thread_data_t *init_job(va_t *adapter_array, adapter_opts_t *opts)
{
    adapter_thread_data_t *job_data = calloc(1, sizeof(*job_data));
    if (!job_data) die("Out of memory\n");

    job_data->adapter_array = copy_adapter_array(adapter_array);
    job_data->adapter = job_data->adapter_array->entries[0];
    job_data->opts = opts;
    job_data->result = -1;
    job_data->next = NULL;

    return job_data;
}

static int processTemplatesThreads(hts_tpool *pool, BAMit_t *bam_in, BAMit_t *bam_out, va_t *adapterArray, adapter_opts_t* opts)
{
    hts_tpool_result *job_result = NULL;
    hts_tpool_process *queue = hts_tpool_process_init(pool, 2 * opts->nthreads, 0);
    adapter_thread_data_t *job_freelist = NULL;
    adapter_thread_data_t *job_data = init_job(adapterArray, opts);
    char qname[257] = { 0 };

    if (!queue) {
        perror("hts_tpool_process_init");
        return -1;
    }

    job_data->record_set = va_init(TEMPLATES_PER_JOB * 2, freeRecord);
    job_data->template_counts = ia_init(TEMPLATES_PER_JOB);
    job_data->nrec = 0;
    job_data->bam_out = bam_out;

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
                int blk = hts_tpool_dispatch2(pool, queue, adapter_job, job_data, 1);
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
                    adapter_thread_data_t *finished_job = hts_tpool_result_data(job_result);
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
                job_data = init_job(adapterArray, opts);
                job_data->record_set = va_init(TEMPLATES_PER_JOB * 2, freeRecord);
                job_data->template_counts = ia_init(TEMPLATES_PER_JOB);
            }
            job_data->nrec = 0;
        }
    }

    if (job_data->template_counts->end > 0) {
        // Deal with left-over items
        if (hts_tpool_dispatch(pool, queue, adapter_job, job_data) < 0) {
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
        adapter_thread_data_t *finished_job = hts_tpool_result_data(job_result);
        output_job_results(bam_out, finished_job);
        finished_job->next = job_freelist;
        job_freelist = finished_job;
        hts_tpool_delete_result(job_result, 0);
    }

    while (job_freelist != NULL) {
        adapter_thread_data_t *next = job_freelist->next;
        job_free(job_freelist);
        job_freelist = next;
    }

    hts_tpool_process_destroy(queue);

    return 0;
}

static int writeMetrics(HashTable *m, adapter_opts_t *opts)
{
    if (!m) return 0;

    HashIter *i;
    HashItem *hi;

    FILE *f = fopen(opts->metrics_name, "w");
    if (!f) { display("Can't open metrics file %s", opts->metrics_name); return -1; }

    i = HashTableIterCreate();
    while ((hi = HashTableIterNext(m, i))) {
        rg_metrics_t *rgm = (rg_metrics_t *)hi->data.p;
        fprintf(f, "RG: %s\n", hi->key);
        fprintf(f, "  Total reads:    %d  %d\n", rgm->total_fwd, rgm->total_rev);
        fprintf(f, "  Adapter reads:  %d  %d\n", rgm->contam_fwd, rgm->contam_rev);
        fprintf(f, "  Histogram Fwd:");
        for (int n = 0; n < rgm->hist_len; n++) {
            if (rgm->hist_fwd[n]) fprintf(f, " %d:%d", n, rgm->hist_fwd[n]);
        }
        fprintf(f, "\n");
        fprintf(f, "  Histogram Rev:");
        for (int n = 0; n < rgm->hist_len; n++) {
            if (rgm->hist_rev[n]) fprintf(f, " %d:%d", n, rgm->hist_rev[n]);
        }
        fprintf(f, "\n");
    }

    HashTableIterDestroy(i);
    fclose(f);
    return 0;
}

/*
 * Main code
 */
static int findAdapters(adapter_opts_t* opts)
{
    int retcode = 1;
    BAMit_t *bam_in = NULL;
    BAMit_t *bam_out = NULL;
    htsThreadPool hts_threads = { NULL, 0 };

    while (1) {
        if (opts->nthreads > 1) {
            hts_threads.pool = hts_tpool_init(opts->nthreads);
            if (!hts_threads.pool) {
                fprintf(stderr, "Couldn't set up thread pool\n");
                break;
            }
        }

        // create metrics hash
        if (opts->metrics_name) {
            metrics = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS);
            if (!metrics) { display("Could not create metrics hash"); break; }
        }

        /*
         * Open input and output BAM files
         */
        bam_in = BAMit_open(opts->input_name, 'r', opts->input_fmt, 0, hts_threads.pool ? &hts_threads : NULL);
        if (!bam_in) break;
        bam_out = BAMit_open(opts->output_name, 'w', opts->output_fmt, opts->compression_level, hts_threads.pool ? &hts_threads : NULL);
        if (!bam_out) break;
        // copy input to output header
        bam_hdr_destroy(bam_out->h); bam_out->h = bam_hdr_dup(bam_in->h);

        // Change header by adding PG lines
        changeHeader(bam_out->h, opts->argv_list);
        if (sam_hdr_write(bam_out->f, bam_out->h) != 0) {
            fprintf(stderr, "Could not write output file header\n");
            break;
        }

        // Read and process each template in the input BAM
        if (opts->nthreads < 2) {
            if (processTemplatesNoThreads(bam_in, bam_out, opts->adapterArray, opts) < 0) break;
        } else {
            if (processTemplatesThreads(hts_threads.pool, bam_in, bam_out, opts->adapterArray, opts) < 0) break;
        }

        if (BAMit_hasnext(bam_in)) break;   // we must has exited the above loop early

        /*
         * And finally.....the metrics
         */
        if (writeMetrics(metrics, opts) != 0) break;
        retcode = 0;

        break;
    }

    // tidy up after us
    BAMit_free(bam_in);
    BAMit_free(bam_out);
    if (hts_threads.pool) hts_tpool_destroy(hts_threads.pool);
    if (metrics) HashTableDestroy(metrics, 1);

    return retcode;
}

/*
 * called from bambi to perform adapter searching
 *
 * Parse the command line arguments, then call the main adapter function
 *
 * returns 0 on success, 1 if there was a problem
 */
int main_adapters(int argc, char *argv[])
{
    int ret = 1;

    adapter_opts_t* opts = parse_args(argc, argv);
    if (opts) {
        ret = findAdapters(opts);
    }
    adapter_free_opts(opts);
    return ret;
}

