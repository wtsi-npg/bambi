/*  read2tags.c -- convert reads into tags.

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

/*
 * position record
 */
typedef struct {
    int record;
    int from;
    int to;
} pos_t;

/*
 * structure to hold options
 */
typedef struct {
    int verbose;
    char *argv_list;
    char compression_level;
    char *in_file;
    char *out_file;
    char *output_fmt;
    char *input_fmt;
    va_t *taglist;
    va_t *qtaglist;
    va_t *poslist;
} opts_t;

/*
 * Release all the options
 */

static void free_opts(opts_t* opts)
{
    if (!opts) return;
    free(opts->argv_list);
    va_free(opts->taglist);
    va_free(opts->qtaglist);
    va_free(opts->poslist);
    free(opts->output_fmt);
    free(opts->input_fmt);
    free(opts->in_file);
    free(opts->out_file);
    free(opts);
}

/*
 * Parse a comma separated list of positions
 * Format is r:s:e,r:s:e,r:s:e,...
 *
 * where r is record number (1 or 2)
 *       s is the start position in the read string
 *       e is the end position in the read string
 * start and end positions are 1 (not zero) based
 */
static void parse_positions(va_t *poslist, char *args)
{
    char *argstr = strdup(args);
    char *save_s;
    char *s = strtok_r(argstr,",",&save_s);
    while (s) {
        pos_t *pos = calloc(1, sizeof(pos_t));
        char *save_p;
        char *p = strtok_r(s,":",&save_p); if (p) pos->record = atoi(p);
        p = strtok_r(NULL,":",&save_p); if (p) pos->from = atoi(p);
        p = strtok_r(NULL,":",&save_p); if (p) pos->to = atoi(p);
        if (!p) {
            // looks like s:e format
            pos->to = pos->from; pos->from = pos->record; pos->record = 0;
        }
        if (pos->record < 0 || pos->record > 2 || pos->from == 0 || pos->to == 0 || pos->from > pos->to) {
            fprintf(stderr,"Invalid pos argument: %s\n", args);
            exit(1);
        }
        va_push(poslist,pos);
        s = strtok_r(NULL,",",&save_s);
    }
    free(argstr);
}

/*
 * display usage information
 */
static void usage(FILE *write_to)
{
    fprintf(write_to,
"Usage: bambi read2tags [options]\n"
"\n"
"Options:\n"
"  -i   --input                 BAM file to read [default: stdin]\n"
"  -o   --output                BAM file to output [default: stdout]\n"
"  -t   --tags                  comma separated list of barcode tags\n"
"  -q   --qtags                 comma separated list of quality  tags\n"
"  -p   --positions             comma separated list of positions\n"
"  -v   --verbose               verbose output\n"
"       --input-fmt             [sam/bam/cram] [default: bam]\n"
"       --output-fmt            [sam/bam/cram] [default: bam]\n"
"       --compression-level     [0..9]\n"
);
}

/*
 * Takes the command line options and turns them into something we can understand
 */
opts_t* read2tags_parse_args(int argc, char *argv[])
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char* optstring = "vi:o:t:q:p:";

    static const struct option lopts[] = {
        { "verbose",            0, 0, 'v' },
        { "input",              1, 0, 'i' },
        { "output",             1, 0, 'o' },
        { "tags",               1, 0, 't' },
        { "qtags",              1, 0, 'q' },
        { "positions",          1, 0, 'p' },
        { "compression-level",  1, 0, 0 },
        { "input-fmt",          1, 0, 0 },
        { "output-fmt",         1, 0, 0 },
        { NULL, 0, NULL, 0 }
    };

    opts_t* opts = calloc(sizeof(opts_t), 1);
    if (!opts) { perror("cannot allocate option parsing memory"); return NULL; }

    opts->argv_list = stringify_argv(argc+1, argv-1);
    if (opts->argv_list[strlen(opts->argv_list)-1] == ' ') opts->argv_list[strlen(opts->argv_list)-1] = 0;

    opts->taglist = va_init(5, free);
    opts->qtaglist = va_init(5, free);
    opts->poslist = va_init(5, free);

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, optstring, lopts, &option_index)) != -1) {
        const char *arg;
        switch (opt) {
        case 'i':   opts->in_file = strdup(optarg);
                    break;
        case 'o':   opts->out_file = strdup(optarg);
                    break;
        case 't':   parse_tags(opts->taglist,optarg);
                    break;
        case 'q':   parse_tags(opts->qtaglist,optarg);
                    break;
        case 'p':   parse_positions(opts->poslist,optarg);
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

    if (!opts->in_file) opts->in_file = strdup("-"); 
    if (!opts->out_file) opts->out_file = strdup("-"); 

    /*
     * And now...validate the inputs
     */
    if (va_isEmpty(opts->taglist)) {
        fprintf(stderr,"You must specify one or more tags\n");
        usage(stderr); return NULL;
    }

    if (va_isEmpty(opts->poslist)) {
        fprintf(stderr,"You must specify one or more positions\n");
        usage(stderr); return NULL;
    }

    if (opts->taglist->end != opts->poslist->end) {
        fprintf(stderr,"You must have the same number of tags and positions\n");
        usage(stderr); return NULL;
    }

    if (opts->taglist->end != opts->qtaglist->end) { 
        fprintf(stderr,"You must have the same number of barcode tags and quality tags\n");
        usage(stderr); return NULL;
    }

    for (int n=0; n < opts->taglist->end; n++) {
        if ( (strlen(opts->taglist->entries[n]) != 2) || (strlen(opts->qtaglist->entries[n]) != 2) ) {
            fprintf(stderr,"Barcode and Quality tags must be two characters\n");
            return NULL;
        }
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
static void sam_hdr_unparse2(SAM_hdr *sh, bam_hdr_t *h)
{
    free(h->text);
    sam_hdr_rebuild(sh);
    h->text = strdup(sam_hdr_str(sh));
    h->l_text = sam_hdr_length(sh);
    sam_hdr_free(sh);
}

/*
 * add PG line to output BAM file
 */
static void addPGLine(BAMit_t *bit, opts_t *opts)
{
    SAM_hdr *sh = sam_hdr_parse_(bit->h->text,bit->h->l_text);

    // add new PG line
    sam_hdr_add_PG(sh, "bambi",
                   "VN", bambi_version(),
                   "CL", opts->argv_list,
                   "DS", "convert reads to tags",
                   NULL, NULL);
    sam_hdr_unparse2(sh,bit->h);
}

/*
 * get and decode the read from a BAM record
 * Why on earth isn't this (and the next few functions) part of htslib?
 * The whole point of an API is that we should not have to know the internal structure...
 */
static char *get_read(bam1_t *rec)
{
    int len = rec->core.l_qseq + 1;
    char *read = calloc(1, kroundup32(len));
    char *seq = bam_get_seq(rec);
    int n;

    for (n=0; n < rec->core.l_qseq; n++) {
        read[n] = "=ACMGRSVTWYHKDBN"[bam_seqi(seq,n)];
    }
    return read;
}

/*
 * get and decode the quality from a BAM record
 */
static char *get_quality(bam1_t *rec)
{
    char *quality = calloc(1, rec->core.l_qseq + 1);
    char *q = bam_get_qual(rec);
    int n;

    for (n=0; n < rec->core.l_qseq; n++) {
        quality[n] = q[n]+33;
    }
    return quality;
}

//
// End of htslib complaints
//

/*
 * create new BAM record from old record with new read and quality strings
 * seq is the uncompressed read string, null terminated
 * qual is the quality string, made printable by adding 33
 * NB this assumes that the new read and qual strings are *not* longer than the originals
 */
static bam1_t *make_new_rec(bam1_t *rec, char *seq, char *qual)
{
    typedef unsigned char uc;
    static const char L[256] = {
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15, 0,15,15,
        15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
        15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
        15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
        15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15
    };
    bam1_t *newrec = bam_dup1(rec);

    newrec->core.l_qseq = strlen(seq);
    int old_len = (rec->core.l_qseq+1)/2 + rec->core.l_qseq;
    int new_len = (newrec->core.l_qseq+1)/2 + newrec->core.l_qseq;
    newrec->l_data = newrec->l_data - old_len + new_len;
;
    char *cp = bam_get_seq(newrec);
    int i;
    for (i = 0; i+1 < newrec->core.l_qseq; i+=2) {
        *cp++ = (L[(uc)seq[i]]<<4) + L[(uc)seq[i+1]];
    }
    if (i < newrec->core.l_qseq)
        *cp++ = L[(uc)seq[i]]<<4;

    if (qual) {
        memcpy(cp, qual, newrec->core.l_qseq);
        for (int n=0; n < newrec->core.l_qseq; n++, *cp++ -= 33);
    } else {
        memset(cp, '\xff', newrec->core.l_qseq);
        cp += newrec->core.l_qseq;
    }

    memcpy(bam_get_aux(newrec), bam_get_aux(rec), bam_get_l_aux(rec));

    return newrec;
}

/*
 * physically remove '\0x01' characters from a string
 */
static void shuffle(char *s)
{
    int len = strlen(s);
    char *p = s;
    while (*p) {
        if (*p == 1) {
            memmove(p,p+1,len);
        } else {
            p++;
        }
        len--;
    }
}

/*
 * Process one record
 */
static bam1_t *process_record(bam1_t *rec, opts_t *opts)
{
    pos_t *pos;
    int recno = -1;
    char *tag_data = calloc(1, rec->core.l_qseq+1);
    char *qtag_data = calloc(1, rec->core.l_qseq+1);

    if (!(rec->core.flag & BAM_FPAIRED)) recno = 0;
    if (rec->core.flag & BAM_FREAD1) recno = 1;
    if (rec->core.flag & BAM_FREAD2) recno = 2;

    char *seq = get_read(rec);
    char *quality = get_quality(rec);

    /*
     * first pass - copy sections of read into tags
     */
    for (int n=0; n < opts->poslist->end; n++) {
        pos = opts->poslist->entries[n];
        if (pos->record == recno) {
            if (pos->from <= rec->core.l_qseq) {
                int from = (pos->from > rec->core.l_qseq) ? rec->core.l_qseq : pos->from;
                int to = (pos->to > rec->core.l_qseq) ? rec->core.l_qseq : pos->to;
                int len = to - from + 1;

                // copy data from read
                memset(tag_data,0,rec->core.l_qseq+1);
                memcpy(tag_data, seq + from - 1, len);
                bam_aux_append(rec, opts->taglist->entries[n], 'Z', len+1, tag_data);

                // copy data from quality
                memset(qtag_data,0,rec->core.l_qseq+1);
                memcpy(qtag_data, quality + from - 1, len);
                bam_aux_append(rec, opts->qtaglist->entries[n], 'Z', len+1, qtag_data);
            }
        }
    }

    /*
     * second pass - mark sections of read as deleted
     */
    for (int n=0; n < opts->poslist->end; n++) {
        pos = opts->poslist->entries[n];
        if (pos->record == recno) {
            if (pos->from <= rec->core.l_qseq) {
                int from = (pos->from > rec->core.l_qseq) ? rec->core.l_qseq : pos->from;
                int to = (pos->to > rec->core.l_qseq) ? rec->core.l_qseq : pos->to;
                int len = to - from + 1;
                memset(seq + from - 1, 1, len);       // mark as deleted
                memset(quality + from - 1, 1, len);   // mark as deleted
            }
        }
    }
    shuffle(seq); shuffle(quality); // physically remove 'marked as deleted' bytes
    bam1_t *newrec = make_new_rec(rec, seq, quality);

    free(tag_data); free(qtag_data); free(quality); free(seq);
    return newrec;
}

/*
 * validate a record before we try to process it
 */
static int invalid_record(bam1_t *rec, int nrec)
{
    if (!(rec->core.flag & BAM_FUNMAP)) {
        fprintf(stderr,"record %d (%s) is aligned. We only handle unaligned records.\n", nrec, bam_get_qname(rec));
        return -1;
    }
    if (rec->core.flag & (BAM_FREVERSE | BAM_FMREVERSE)) {
        fprintf(stderr,"record %d (%s) is reversed. We can't handle that.\n", nrec, bam_get_qname(rec));
        return -1;
    }
    if (rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
        fprintf(stderr,"record %d (%s) is secondary or supplementary. We can't handle that.\n", nrec, bam_get_qname(rec));
        return -1;
    }
    return 0;
}

/*
 * return the length of some aux data
 */
static int aux_type2size(char *s)
{
    switch (*s) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return strlen(s+1) + 1;
    default:
        return 0;
    }
}

/*
 * merge two read records, where one of them has no reads
 */
static bam1_t *merge_records(bam1_t *r1, bam1_t *r2)
{
    bam1_t *src, *dst;

    if (r1->core.l_qseq==0 && r2->core.l_qseq==0) {
        fprintf(stderr,"Both records are empty (%s) - aborting\n", bam_get_qname(r1));
        exit(1);
    }

    // start with the non-empty record...
    if (r1->core.l_qseq) { dst = bam_dup1(r1); src = r2; }
    else                 { dst = bam_dup1(r2); src = r1; }

    // copy aux tags from src to dst
    uint8_t *s;
    s = bam_get_aux(src);
    while ((s - bam_get_aux(src)) < bam_get_l_aux(src)) {
        char tag[3];
        tag[0] = *s++;
        tag[1] = *s++;
        tag[2] = 0;
        char type = *s;
        int len = aux_type2size(s++);
        bam_aux_append(dst,tag,type,len,s);
        s += len;
    }

    // turn of paired flags to make single read
    dst->core.flag &= ~(BAM_FREAD1 | BAM_FREAD2 | BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FMUNMAP);
    return dst;
}

/*
 * write BAM record, and check for failure
 */
static int write_record(BAMit_t *bam, bam1_t *rec)
{
    int r = sam_write1(bam->f, bam->h, rec);
    if (r < 0) {
        fprintf(stderr,"sam_write1() failed\n");
        return -1;
    }
    return 0;
}

        
/*
 * Main code
 *
 * Open and process the BAM file
 */
int process(opts_t* opts)
{
    int retcode = 1;
    int nrec = 0;
    int r;

    BAMit_t *bam_in = BAMit_open(opts->in_file, 'r', opts->input_fmt, 0);
    BAMit_t *bam_out = BAMit_open(opts->out_file, 'w', opts->output_fmt, opts->compression_level);

    // copy input to output header
    bam_hdr_destroy(bam_out->h); bam_out->h = bam_hdr_dup(bam_in->h);
    addPGLine(bam_out, opts);               // Add a new PG line
    r = sam_hdr_write(bam_out->f, bam_out->h);  // and write header
    if (r<0) {
        fprintf(stderr,"Can't write header\n");
        return -1;
    }

    while (BAMit_hasnext(bam_in)) {
        bam1_t *rec = BAMit_next(bam_in);
        if (invalid_record(rec,++nrec)) return -1;
        bam1_t *newrec = process_record(rec,opts);

        bam1_t *rec2 = BAMit_peek(bam_in);
        if (rec2 && strcmp(bam_get_qname(rec), bam_get_qname(rec2)) == 0) {
            rec2 = BAMit_next(bam_in);
            if (invalid_record(rec2,++nrec)) return -1;
            bam1_t *newrec2 = process_record(rec2,opts);
            if ((newrec->core.l_qseq == 0) || (newrec2->core.l_qseq == 0)) {
                bam1_t *merged_rec = merge_records(newrec, newrec2);
                if (write_record(bam_out, merged_rec)) return -1;
                bam_destroy1(merged_rec);
            } else {
                if (write_record(bam_out, newrec)) return -1;
                if (write_record(bam_out, newrec2)) return -1;
            }
            bam_destroy1(newrec2);
        } else {
            if (write_record(bam_out,newrec)) return -1;
        }
        bam_destroy1(newrec);
    }

    // tidy up after us
    BAMit_free(bam_in);
    BAMit_free(bam_out);
    return retcode;
}

/*
 * called from bambi to convert reads to tags
 *
 * Parse the command line arguments, then call the main process() function
 *
 * returns 0 on success, 1 if there was a problem
 */
int main_read2tags(int argc, char *argv[])
{
    int ret = 1;
    opts_t* opts = read2tags_parse_args(argc, argv);
    if (opts) ret = process(opts);
    free_opts(opts);
    return ret;
}
