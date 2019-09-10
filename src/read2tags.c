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
#include <math.h>

#include "array.h"
#include "bamit.h"
#include "parse.h"

#define DEFAULT_KEEP_TAGS "BC,QT,RG"
#define DEFAULT_DISCARD_TAGS "as,af,aa,a3,ah"


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
    va_t *keep_tags;
    va_t *discard_tags;
    bool merge;
    bool replace;
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
    va_free(opts->keep_tags);
    va_free(opts->discard_tags);
    free(opts->output_fmt);
    free(opts->input_fmt);
    free(opts->in_file);
    free(opts->out_file);
    free(opts);
}

/*
 * Compare two aux tags for equality.
 * NB: this needs to move into htslib....
 */
static int bam_aux_cmp(const uint8_t *s, const uint8_t *d)
{
    if (*s == *d) {
        // same type
        char type = *s;
        s++; d++;   // skip over type
        switch (type) {
            case 'c': if ((int32_t)*(int8_t*)s == (int32_t)*(int8_t*)d) return 0;
                      break;
            case 'C': if ((int32_t)*(uint8_t*)s == (int32_t)*(uint8_t*)d) return 0;
                      break;
            case 's': if ((int32_t)*(int16_t*)s == (int32_t)*(int16_t*)d) return 0;
                      break;
            case 'S': if ((int32_t)*(uint16_t*)s == (int32_t)*(uint16_t*)d) return 0;
                      break;
            case 'i':
            case 'I': if (*(int32_t*)s == *(int32_t*)d) return 0;
                      break;
            case 'A': if (*(char*)s == *(char*)d) return 0;
                      break;
            case 'd': if (fabs(*(double*)s - *(double*)d) < .0001) return 0;
                      break;
            case 'f': if (fabsf(*(float*)s - *(float*)d) < .0001) return 0;
                      break;
            case 'Z':
            case 'H': if (strcmp((char*)s, (char *)d) == 0) return 0;
                      break;
            default: fprintf(stderr,"bam_aux_cmp(): unrecognized type [%c]\n", type);
        }
    }
    return 1;
}

/*
 * Parse a comma separated list of positions
 * Format is r:s:e,r:s:e,r:s:e,...
 *
 * where r is record number (0, 1 or 2 and is optional)
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
"  -k   --keep-tags             comma separated list of tags to keep when merging records\n"
"                               [default: " DEFAULT_KEEP_TAGS "]\n"
"  -d   --discard-tags          comma separated list of tags to discard when merging records\n"
"                               [default: " DEFAULT_DISCARD_TAGS "]\n"
"  -p   --positions             comma separated list of positions\n"
"  -m   --merge                 merge duplicate tags\n"
"  -r   --replace               replace duplicate tags\n"
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

    const char* optstring = "vi:o:t:q:p:k:d:mr";

    static const struct option lopts[] = {
        { "verbose",            0, 0, 'v' },
        { "input",              1, 0, 'i' },
        { "output",             1, 0, 'o' },
        { "tags",               1, 0, 't' },
        { "qtags",              1, 0, 'q' },
        { "positions",          1, 0, 'p' },
        { "keep-tags",          1, 0, 'k' },
        { "discard-tags",       1, 0, 'd' },
        { "merge",              0, 0, 'm' },
        { "replace",            0, 0, 'r' },
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
    opts->keep_tags = va_init(10, free);
    opts->discard_tags = va_init(10, free);

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
        case 'k':   parse_tags(opts->keep_tags,optarg);
                    break;
        case 'd':   parse_tags(opts->discard_tags,optarg);
                    break;
        case 'v':   opts->verbose++;
                    break;
        case 'm':   opts->merge = true;
                    break;
        case 'r':   opts->replace = true;
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

    if (va_isEmpty(opts->keep_tags)) parse_tags(opts->keep_tags,DEFAULT_KEEP_TAGS);
    if (va_isEmpty(opts->discard_tags)) parse_tags(opts->discard_tags,DEFAULT_DISCARD_TAGS);
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

    if (opts->merge && opts->replace) {
        fprintf(stderr,"You can't have --merge AND --replace. Choose one.\n");
        usage(stderr); return NULL;
    }

    return opts;
}

/*
 * add PG line to output BAM file
 */
static void addPGLine(BAMit_t *bit, opts_t *opts)
{
    bam_hdr_t *hdr = bit->h;

    // add new PG line
    sam_hdr_add_pg(hdr, "bambi",
                   "VN", bambi_version(),
                   "CL", opts->argv_list,
                   "DS", "convert reads to tags",
                   NULL, NULL);
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
    uint8_t *seq = bam_get_seq(rec);
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
    uint8_t *q = bam_get_qual(rec);
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
    uint8_t *cp = bam_get_seq(newrec);
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
 * add a new tag to our taglist, or append to existing tag
 */
static void add_or_update(va_t *va, char *tag, char *data)
{
    int n;
    for (n=0; n < va->end; n++) {
        if (strncmp(tag,va->entries[n],2) == 0) break;
    }

    if (n == va->end) {
        // add new tag
        char *e = calloc(1, strlen(tag) + 1 + strlen(data) + 1);
        strcpy(e, tag);
        strcat(e, ":");
        strcat(e, data);
        va_push(va,e);
    } else {
        // update existing tag
        va->entries[n] = realloc(va->entries[n], strlen(va->entries[n]) + strlen(data) + 1);
        strcat(va->entries[n],data);
    }
    return;
}

/*
 * add new tag to the record, or replace or update if it already exists
 */
static void add_tag(bam1_t *rec, char *tag, char *data, opts_t *opts)
{
    uint8_t *s = bam_aux_get(rec,tag);
    if (s) { // tag already exists
        if (opts->replace) {
            bam_aux_del(rec,s);
            bam_aux_append(rec, tag, 'Z', strlen(data)+1, (uint8_t *) data);
        }
        if (opts->merge) {
            if (*s != 'Z' && *s != 'H') { 
                fprintf(stderr,"Trying to merge tag [%s] which is type [%c]\n", tag, *s);
                exit(1);
            }
            char *old_data = bam_aux2Z(s);
            size_t old_len = strlen(old_data);
            size_t data_len = strlen(data);
            uint8_t *new_data = malloc(old_len+data_len+1);
            if (!new_data) die("Out of memory");
            memcpy(new_data, old_data, old_len);
            memcpy(new_data + old_len, data, data_len + 1);
            bam_aux_del(rec,s);
            bam_aux_append(rec, tag, 'Z', old_len+data_len+1, new_data);
            free(new_data);
        }
        if (!opts->replace && !opts->merge) {
            fprintf(stderr,"Found duplicate tag [%s] and no --replace or --merge option\n", tag);
            exit(1);
        }
    } else {
        // add new tag
        bam_aux_append(rec, tag, 'Z', strlen(data)+1, (uint8_t *) data);
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
    va_t *new_tags = va_init(10,free);
    va_t *new_qtags = va_init(10,free);

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
                add_or_update(new_tags, opts->taglist->entries[n], tag_data);

                // copy data from quality
                memset(qtag_data,0,rec->core.l_qseq+1);
                memcpy(qtag_data, quality + from - 1, len);
                add_or_update(new_qtags, opts->qtaglist->entries[n], qtag_data);

            }
        }
    }

    // add new tags
    for (int n=0; n < new_tags->end; n++) {
        char *tag = new_tags->entries[n]; tag[2] = 0;
        char *data = tag+3;
        add_tag(rec, tag, data, opts);
    }

    // add new quality tags
    for (int n=0; n < new_qtags->end; n++) {
        char *tag = new_qtags->entries[n]; tag[2] = 0;
        char *data = tag+3;
        add_tag(rec, tag, data, opts);
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
    va_free(new_tags); va_free(new_qtags);
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
static int aux_type2size(uint8_t *s)
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
        return strlen((char *) s+1) + 1;
    default:
        return 0;
    }
}

/*
 * merge two read records, where one of them has no reads
 */
static bam1_t *merge_records(bam1_t *r1, bam1_t *r2, opts_t *opts)
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

        // is the tag in our 'keep' list? (or taglist or qtaglist?)
        if ( (va_contains(opts->keep_tags, tag) != -1) ||
             (va_contains(opts->taglist, tag) != -1) ||
             (va_contains(opts->qtaglist, tag) != -1) ) {

            if (bam_aux_get(dst, tag) == 0) {   // not already there
                bam_aux_append(dst,tag,type,len,s);
            } else {
                if (opts->merge) {
                    // merge with existing tag
                    if (type == 'Z' || type == 'H') {
                        uint8_t *t = bam_aux_get(dst,tag);
                        size_t t_len = strlen((char *) t + 1);
                        size_t s_len = strlen((char *) s);
                        uint8_t *data = malloc(t_len + s_len + 1);
                        if (!data) die("Out of memory");
                        memcpy(data, t + 1, t_len);
                        memcpy(data + t_len, s, s_len + 1);
                        bam_aux_del(dst,t);
                        bam_aux_append(dst,tag,type,t_len + s_len + 1,data);
                        free(data);
                    }
                }
                if (opts->replace) {
                    // replace existing tag
                    uint8_t *t = bam_aux_get(dst,tag);
                    bam_aux_del(dst,t);
                    bam_aux_append(dst,tag,type,len,s);
                }
                if (!opts->merge && !opts->replace) {
                    uint8_t *t = bam_aux_get(dst,tag);
                    if (bam_aux_cmp(s-1,t)) {
                        fprintf(stderr,"Tag [%s] already exists and is not the same value\n", tag);
                        exit(1);
                    }
                }
            }
        } else {
            if (va_contains(opts->discard_tags, tag) == -1) {
                fprintf(stderr,"Tag %s is in neither keep nor discard list\n", tag);
                exit(1);
            }
        }

        s += len;
    }

    // turn off paired flags to make single read
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
    int retcode = 0;
    int nrec = 0;
    int r;

    BAMit_t *bam_in = BAMit_open(opts->in_file, 'r', opts->input_fmt, 0, NULL);
    BAMit_t *bam_out = BAMit_open(opts->out_file, 'w', opts->output_fmt, opts->compression_level, NULL);

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
                bam1_t *merged_rec = merge_records(newrec, newrec2, opts);
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
