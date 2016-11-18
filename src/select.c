/*  select.c -- select alignments into separate BAM files.

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

void closeSamFile(void *f) { hts_close((samFile *)f); }
void freeBamHdr(void *h) { if (h) bam_hdr_destroy((bam_hdr_t *)h); h=NULL; }
void freeRecord(void *r) { bam_destroy1((bam1_t *)r); }
void freeRecordSet(void *s) { va_free((va_t *)s); }
void freeChimera(void *s) { ia_free((ia_t *)s); }

bool verbose;

/*
 * structure to hold options
 */
typedef struct {
    int verbose;
    char *argv_list;
    char compression_level;
    va_t *in_file;
    va_t *out_file;
    char *unaligned_filename;
    char *metrics_filename;
    char *output_fmt;
    char *input_fmt;
} opts_t;

/****
void freeSQList(void *sq)
{
    va_free(sq);
}
*/

/*
 * structure to hold metrics
 */
typedef struct {
    int nAlignments;
    int nReads;
    int nUnaligned;
    ia_t *nReadsPerRef;
    va_t *nChimericReads;
    ia_t *nAlignedForward;
    ia_t *nAlignedReverse;
} metrics_t;

metrics_t *metrics_init(int sz)
{
    metrics_t *m = calloc(sizeof(metrics_t), 1);
    m->nReads = 0;
    m->nUnaligned = 0;
    m->nReadsPerRef = ia_init(sz);
    m->nChimericReads = va_init(sz, freeChimera);
    m->nAlignedForward = ia_init(sz+1);
    m->nAlignedReverse = ia_init(sz+1);

    m->nAlignments = sz;
    int n;
    for (n=0; n < sz; n++) {
        ia_push(m->nReadsPerRef,0);
        ia_push(m->nAlignedForward,0);
        ia_push(m->nAlignedReverse,0);
        ia_t *ia = ia_init(sz);
        int i;
        for (i=0; i<sz; i++) { ia_push(ia,0); }
        va_push(m->nChimericReads,ia);
    }
    ia_push(m->nAlignedForward,0);
    ia_push(m->nAlignedReverse,0);
    return m;
}

void metrics_free(metrics_t *m)
{
    if (!m) return;
    ia_free(m->nReadsPerRef);
    va_free(m->nChimericReads);
    ia_free(m->nAlignedForward);
    ia_free(m->nAlignedReverse);
    free(m);
}


/*
 * stolen from samtools/padding.c
 */
static void replace_cigar(bam1_t *b, int n, uint32_t *cigar)
{
    if (n != b->core.n_cigar) {
        int o = b->core.l_qname + b->core.n_cigar * 4;
        if (b->l_data + (n - b->core.n_cigar) * 4 > b->m_data) {
            b->m_data = b->l_data + (n - b->core.n_cigar) * 4;
            kroundup32(b->m_data);
            b->data = (uint8_t*)realloc(b->data, b->m_data);
        }
        memmove(b->data + b->core.l_qname + n * 4, b->data + o, b->l_data - o);
        memcpy(b->data + b->core.l_qname, cigar, n * 4);
        b->l_data += (n - b->core.n_cigar) * 4;
        b->core.n_cigar = n;
    } else memcpy(b->data + b->core.l_qname, cigar, n * 4);
}


/*
 * Release all the options
 */

static void free_opts(opts_t* opts)
{
    if (!opts) return;
    free(opts->argv_list);
    va_free(opts->in_file);
    va_free(opts->out_file);
    free(opts->unaligned_filename);
    free(opts->metrics_filename);
    free(opts->output_fmt);
    free(opts->input_fmt);
    free(opts);
}

/*
 * display usage information
 */
static void usage(FILE *write_to)
{
    fprintf(write_to,
"Usage: bambi select [options]\n"
"\n"
"Options:\n"
"  -i   --input                 comma separated list of BAM files to read\n"
"  -o   --output                comma separated list of BAM files to output\n"
"  -n                           BAM file to write unaligned reads to [optional]\n"
"  -m                           file to write metrics to [optional]\n"
"  -v   --verbose               verbose output\n"
"       --input-fmt             [sam/bam/cram] [default: bam]\n"
"       --output-fmt            [sam/bam/cram] [default: bam]\n"
"       --compression-level     [0..9]\n"
);
}

/*
 * Takes the command line options and turns them into something we can understand
 */
opts_t* select_parse_args(int argc, char *argv[])
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char* optstring = "vi:o:n:m:";

    static const struct option lopts[] = {
        { "verbose",            0, 0, 'v' },
        { "input",              1, 0, 'i' },
        { "output",             1, 0, 'o' },
        { "compression-level",  1, 0, 0 },
        { "input-fmt",          1, 0, 0 },
        { "output-fmt",         1, 0, 0 },
        { NULL, 0, NULL, 0 }
    };

    opts_t* opts = calloc(sizeof(opts_t), 1);
    if (!opts) { perror("cannot allocate option parsing memory"); return NULL; }

    opts->argv_list = stringify_argv(argc+1, argv-1);
    if (opts->argv_list[strlen(opts->argv_list)-1] == ' ') opts->argv_list[strlen(opts->argv_list)-1] = 0;

    opts->in_file = va_init(5, free);
    opts->out_file = va_init(5, free);

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, optstring, lopts, &option_index)) != -1) {
        const char *arg;
        switch (opt) {
        case 'i':   parse_tags(opts->in_file, optarg);
                    break;
        case 'o':   parse_tags(opts->out_file, optarg);
                    break;
        case 'n':   opts->unaligned_filename = strdup(optarg);
                    break;
        case 'm':   opts->metrics_filename = strdup(optarg);
                    break;
        case 'v':   opts->verbose++;
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

    if (va_isEmpty(opts->in_file)) {
        fprintf(stderr,"You must specify one or more input files\n");
        usage(stderr); return NULL;
    }

    if (va_isEmpty(opts->out_file)) {
        fprintf(stderr,"You must specify one or more output files\n");
        usage(stderr); return NULL;
    }

    if (opts->in_file->end != opts->out_file->end) {
        fprintf(stderr,"You must have the same number of input and output files\n");
        usage(stderr); return NULL;
    }

    if (opts->compression_level && !isdigit(opts->compression_level)) {
        fprintf(stderr, "compression-level must be a digit in the range [0..9], not '%c'\n", opts->compression_level);
        usage(stderr); return NULL;
    }

    return opts;
}

/*
 * convert SAM_hdr to bam_hdr
 */
void sam_hdr_unparse2(SAM_hdr *sh, bam_hdr_t *h)
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
void addHeaderLines(BAMit_t *bit, opts_t *opts, bool u)
{
    SAM_hdr *sh = sam_hdr_parse_(bit->h->text,bit->h->l_text);

    if (u) {
        // for the unaligned BAM file, we need to remove the *SQ lines
        sh = sam_hdr_del(sh, "SQ", NULL, NULL);
    }

    // specify sort order
    sh->sort_order = ORDER_UNSORTED;

    // add new PG line
    sam_hdr_add_PG(sh, "bambi",
                   "VN", bambi_version(),
                   "CL", opts->argv_list,
                   "DS", "Split alignments into different files",
                   NULL, NULL);
    sam_hdr_unparse2(sh,bit->h);
}

/*
 * Open ALL the files
 */
int openSamFiles(va_t *in_file, va_t *out_file, va_t *in_bit, va_t *out_bit, opts_t *opts)
{
    int n;
    samFile *f;
    bam_hdr_t *h;
    BAMit_t *bit;

    for (n=0; n < in_file->end; n++) {
        // open input file
        f = openSamFile(in_file->entries[n], opts->input_fmt, opts->compression_level, 'r');
        h = sam_hdr_read(f);
        bit = BAMit_init(f,h);
        va_push(in_bit,bit);

        // open output file
        f = openSamFile(out_file->entries[n], opts->output_fmt, opts->compression_level, 'w');
        bit = BAMit_init(f,bam_hdr_dup(h));
        va_push(out_bit,bit);

        // add PG line
        addHeaderLines(bit,opts,false);

        // write header
        if (sam_hdr_write(bit->f,bit->h)) { fprintf(stderr,"Failed to write header\n"); exit(1); }
    }
    return 0;
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
 * Look through the list of record sets to find the first record set containing an aligned read
 */
static int firstAlignedIndex(va_t *recordSetList)
{
    int setnum, recnum;

    // for each record set
    for (setnum=0; setnum < recordSetList->end; setnum++) {
        va_t *recordSet = recordSetList->entries[setnum];
        // look for an aligned read in the set
        for (recnum=0; recnum < recordSet->end; recnum++) {
            bam1_t *rec = recordSet->entries[recnum];
            if (!(rec->core.flag & BAM_FUNMAP)) {
                return setnum;      // this set contains an aligned read
            }
        }
    }
    return -1;  // no alignments found
}

/*
 * write a record set to an output BAM file
 */
void writeRecordSet(BAMit_t *bit, va_t *recordSet)
{
    int n,r;

    for (n=0; n < recordSet->end; n++) {
        bam1_t *rec = recordSet->entries[n];
        if (rec->core.flag & BAM_FUNMAP) {
            // unmapped! So reset values
            rec->core.tid=-1;
            rec->core.mtid=-1;
            rec->core.pos=-1;
            rec->core.mpos=-1;
            rec->core.qual=0;
            rec->core.isize=0;
            replace_cigar(rec,0,(uint32_t *)"*");
            uint8_t *s = bam_aux_get(rec,"MD");
            if (s) bam_aux_del(rec,s);
        }
        r = sam_write1(bit->f, bit->h, rec);
        if (r <= 0) {
            fprintf(stderr, "Problem writing record %d : %d\n", n,r);
            exit(1);
        }
    }
}

/*
 * write references to metrics file
 */
void writeReferences(FILE *f, va_t *in_bit)
{
    int n, nSQ;

    fprintf(f,"\"refList\":[");
    for (n=0; n < in_bit->end; n++) {
        BAMit_t *bit = in_bit->entries[n];
        SAM_hdr *sh = sam_hdr_parse_(bit->h->text,bit->h->l_text);
        if (n) fprintf(f,",");
        fprintf(f,"[");
        for (nSQ=0; nSQ<sh->nref; nSQ++) {
            SAM_SQ sq = sh->ref[nSQ];
            if (nSQ) fprintf(f,",");
            fprintf(f,"{");
            SAM_hdr_tag *tag = sq.tag;
            char *ur=NULL, *ln=NULL, *sp=NULL, *as=NULL, *sn=NULL;
            while (tag) {
                if (strlen(tag->str)>3) {
                    if (strncasecmp(tag->str,"UR:",3) == 0) ur = strdup(tag->str+3);
                    if (strncasecmp(tag->str,"LN:",3) == 0) ln = strdup(tag->str+3);
                    if (strncasecmp(tag->str,"SP:",3) == 0) sp = strdup(tag->str+3);
                    if (strncasecmp(tag->str,"AS:",3) == 0) as = strdup(tag->str+3);
                    if (strncasecmp(tag->str,"SN:",3) == 0) sn = strdup(tag->str+3);
                }
                tag = tag->next;
            }
            if (ur) fprintf(f,"\"ur\":\"%s\",", ur);   else fprintf(f,"\"ur\":null,");
            if (ln) fprintf(f,"\"ln\":%d,", atoi(ln)); else fprintf(f,"\"ln\":null,");
            if (sp) fprintf(f,"\"sp\":\"%s\",", sp);   else fprintf(f,"\"sp\":null,");
            if (as) fprintf(f,"\"as\":\"%s\",", as);   else fprintf(f,"\"as\":null,");
            if (sn) fprintf(f,"\"sn\":\"%s\"", sn);    else fprintf(f,"\"sn\":null");
            fprintf(f,"}");
            free(ur); free(ln); free(sp); free(as); free(sn);
        }
        fprintf(f,"]");
        sam_hdr_free(sh);
    }
    fprintf(f,"],");
}

/*
 * Write the metrics file
 */
static void writeMetrics(va_t *in_bit, metrics_t *metrics, opts_t *opts)
{
    char *s;
    // Open the metrics file
    FILE *f = fopen(opts->metrics_filename, "w");
    if (!f) {
        fprintf(stderr,"Can't create metrics file %s\n", opts->metrics_filename);
        exit(1);
    }

    fprintf(f,"{");

    // references
    writeReferences(f,in_bit);

    // Chimeric read counts
    int n;
    fprintf(f,"\"chimericReadsCount\":[");
    for (n=0; n < metrics->nChimericReads->end; n++) {
        ia_t *ia = metrics->nChimericReads->entries[n];
        s = ia_join(ia,",");
        if (n) fprintf(f,",");
        fprintf(f,"[%s]",s);
        free(s);
    }
    fprintf(f,"],");

    s = ia_join(metrics->nAlignedReverse,",");
    fprintf(f,"\"readsCountByAlignedNumReverse\":[%s],",s); free(s);
    s = ia_join(metrics->nAlignedForward,",");
    fprintf(f,"\"readsCountByAlignedNumForward\":[%s],",s); free(s);
    s = ia_join(metrics->nReadsPerRef,",");
    fprintf(f,"\"readsCountPerRef\":[%s],", s); free(s);
    fprintf(f,"\"readsCountUnaligned\":%d,", metrics->nUnaligned);
    fprintf(f,"\"totalReads\":%d,", metrics->nReads);
    fprintf(f,"\"numberAlignments\":%d,", metrics->nAlignments);
    fprintf(f,"\"programVersion\":\"%s\",", bambi_version());
    fprintf(f,"\"programCommand\":\"%s\",", opts->argv_list);
    fprintf(f,"\"programName\":\"bambi\"");
    fprintf(f,"}");
    fclose(f);
}

ia_t *checkAlignmentsByRef(va_t *recordSetList, bool result)
{
    ia_t *ia = ia_init(5);
    int i,j;

    for (i=0; i < recordSetList->end; i++) {
        int found = 0;
        va_t *recordSet = recordSetList->entries[i];
        for (j=0; j < recordSet->end; j++) {
            bam1_t *rec = recordSet->entries[j];
            if ( (rec->core.flag & BAM_FPAIRED) && ( (bool)(rec->core.flag & BAM_FREAD2) == result )) {
                if ( !(rec->core.flag & BAM_FUNMAP) ) {
                    found = 1;
                }
            }
        }
        ia_push(ia,found);
    }
    return ia;
}

static int indexAlignment(ia_t *ia)
{
    int n;
    for (n=0; n < ia->end; n++) {
        if (ia->entries[n]) return n;
    }
    return -1;
}

void checkNextReadsForChimera(va_t *recordSetList, metrics_t *metrics)
{
    ia_t *alignmentByRef = checkAlignmentsByRef(recordSetList,false);
    ia_t *alignmentByRefPaired = checkAlignmentsByRef(recordSetList,true);

    int sumAlignments = ia_sum(alignmentByRef);
    int sumAlignmentsPaired = ia_sum(alignmentByRefPaired);

    if (sumAlignments == 1 && sumAlignmentsPaired == 1) {
        int indexRef = indexAlignment(alignmentByRef);
        int indexRefPaired = indexAlignment(alignmentByRefPaired);

        if (indexRef != -1 && indexRefPaired != -1) {
           // getChimericReadsCount()[indexRef][indexRefPaired]++;
            ia_t *ia = metrics->nChimericReads->entries[indexRef];
            ia->entries[indexRefPaired]++;
        }

        assert(indexRef == indexRefPaired);
    }

    metrics->nAlignedForward->entries[sumAlignments]++;
    metrics->nAlignedReverse->entries[sumAlignmentsPaired]++;

    ia_free(alignmentByRefPaired);
    ia_free(alignmentByRef);
}

/*
 * Read the array of input files, write to the array of output files
 * (and the optional unaligned_bam file)
 */
static int processFiles(va_t *in_bit, va_t *out_bit, BAMit_t *unaligned_bam, opts_t *opts)
{
    BAMit_t *outBam;
    BAMit_t *bit;
    metrics_t *metrics = metrics_init(in_bit->end);
    int n;


    BAMit_t *firstBit = in_bit->entries[0];
    while (BAMit_hasnext(firstBit)) {
        metrics->nReads++;
        // for each record set
        va_t *recordSetList = va_init(5,freeRecordSet);
        bam1_t *rec = BAMit_peek(firstBit);
        char *qname = strdup(bam_get_qname(rec));
        for (n=0; n < in_bit->end; n++) {
            // for each BAM file
            va_t *recordSet = read_record_set(in_bit->entries[n],qname);
            va_push(recordSetList,recordSet);
        }
        free(qname);

        checkNextReadsForChimera(recordSetList, metrics);

        n = firstAlignedIndex(recordSetList);
        if (n == -1) {
            // unaligned
            outBam = unaligned_bam ? unaligned_bam : out_bit->entries[out_bit->end-1];
            n = recordSetList->end-1;
            metrics->nUnaligned++;
        } else {
            outBam = out_bit->entries[n];
            metrics->nReadsPerRef->entries[n]++;
        }
        writeRecordSet(outBam,recordSetList->entries[n]);
        va_free(recordSetList);
    }

    if (opts->metrics_filename) writeMetrics(in_bit, metrics, opts);
    metrics_free(metrics);

    return 0;
}

/*
 * Main code
 *
 * Open all the BAM files, then process them
 */
int aln_select(opts_t* opts)
{
    int retcode = 1;

    BAMit_t *unaligned_bam = NULL;
    va_t *in_bit = va_init(5,BAMit_free);
    va_t *out_bit = va_init(5,BAMit_free);

    openSamFiles(opts->in_file, opts->out_file, in_bit, out_bit, opts);

    if (opts->unaligned_filename) {
        samFile *f = openSamFile(opts->unaligned_filename, opts->output_fmt, opts->compression_level, 'w');
        if (!f) {
            fprintf(stderr,"Can't open unaligned output file (%s)\n", opts->unaligned_filename);
            return 1;
        }
        BAMit_t *bit = in_bit->entries[in_bit->end-1];
        bam_hdr_t *h = bam_hdr_dup(bit->h);
        unaligned_bam = BAMit_init(f,h);
        addHeaderLines(unaligned_bam,opts,true);
        int nt = unaligned_bam->h->n_targets;  // Aaaargh!
        unaligned_bam->h->n_targets=0;  // Aaaargh!
        if (sam_hdr_write(f,unaligned_bam->h)) { fprintf(stderr,"Failed to write header\n"); exit(1); }
        unaligned_bam->h->n_targets=nt;  // Aaaargh!
    }

    if (in_bit && out_bit) retcode = processFiles(in_bit, out_bit, unaligned_bam, opts);

    // tidy up after us
    va_free(in_bit);
    va_free(out_bit);
    BAMit_free(unaligned_bam);

    return retcode;
}

/*
 * called from bambi to perform selection by alignment
 *
 * Parse the command line arguments, then call the main select() function
 *
 * returns 0 on success, 1 if there was a problem
 */
int main_select(int argc, char *argv[])
{
    int ret = 1;
    opts_t* opts = select_parse_args(argc, argv);
    verbose = opts->verbose;
    if (opts) ret = aln_select(opts);
    free_opts(opts);
    return ret;
}
