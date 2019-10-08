/*  bamit.c -- iterator for reading BAM files one record at a time.

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

#include <stdio.h>

#include "bambi_utils.h"
#include "bamit.h"

BAMit_t *BAMit_init(samFile *f, sam_hdr_t *h)
{
    int r;
    BAMit_t *bit = calloc(1,sizeof(BAMit_t));
    bit->f = f;
    bit->h = h;
    bit->rec = bam_init1();
    bit->nextRec = bam_init1();
    if (f->is_write == 0) {
        r = sam_read1(bit->f, bit->h, bit->nextRec);
        if (r<0) { bam_destroy1(bit->nextRec); bit->nextRec = NULL; }
    }
    return bit;
}

/*
 * Open a BAM file
 * arguments are: char *fname               filename to open
 *                char mode                 'r' or 'w'
 *                char *fmt                 format [bam,sam,cram]
 *                char compression level    [0..9]
 */
BAMit_t *BAMit_open(char *fname, char mode, char *fmt, char compression_level,
                    htsThreadPool *thread_pool)
{
    samFile *f = NULL;
    sam_hdr_t *h = NULL;
    htsFormat *format = NULL;
    char m[5];

    if (fmt) {
        format = calloc(1,sizeof(htsFormat));
        if (hts_parse_format(format, fmt) < 0) {
            fprintf(stderr,"Unknown input format: %s\n", fmt);
            exit(1);
        }
    }
    sam_open_mode(m+1, fname, NULL);    // set type (sam/bam/cram) from filename
    m[0] = mode;
    m[2] = compression_level;
    f = hts_open_format(fname, m, format);
    free(format);
    if (!f) {
        fprintf(stderr,"Could not open file (%s)\n", fname);
        exit(1);
    }

    if (thread_pool && hts_set_thread_pool(f, thread_pool) < 0) {
        fprintf(stderr, "Couldn't set thread pool on %s\n", fname);
        exit(1);
    }

    if (mode == 'r') h = sam_hdr_read(f);
    else             h = bam_hdr_init();

    return BAMit_init(f,h);
}

void BAMit_free(void *ptr)
{
    BAMit_t *bit = (BAMit_t *)ptr;
    if (!bit) return;
    if (bit->f) hts_close(bit->f);
    if (bit->h) bam_hdr_destroy(bit->h);
    if (bit->rec) bam_destroy1(bit->rec);
    if (bit->nextRec) bam_destroy1(bit->nextRec);
    free(bit);
}

bam1_t *BAMit_next(BAMit_t *bit)
{
    if (!bit->nextRec) return NULL;
    if (!bam_copy1(bit->rec,bit->nextRec)) die("bam_copy1() failed in BAMit_next()");
    int r = sam_read1(bit->f, bit->h, bit->nextRec);
    if (r<0) { bam_destroy1(bit->nextRec); bit->nextRec = NULL; }
    return bit->rec;
}

bam1_t *BAMit_peek(BAMit_t *bit)
{
    return bit->nextRec;
}

bool BAMit_hasnext(BAMit_t *bit)
{
    return (bit->nextRec != NULL);
}




