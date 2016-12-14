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

#include "bamit.h"

BAMit_t *BAMit_init(samFile *f, bam_hdr_t *h)
{
    int r;
    BAMit_t *bit = calloc(1,sizeof(BAMit_t));
    bit->f = f;
    bit->h = h;
    bit->rec = bam_init1();
    bit->nextRec = bam_init1();
    if (f->is_write == 0) r = sam_read1(bit->f, bit->h, bit->nextRec);
    return bit;
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
    bam_copy1(bit->rec,bit->nextRec);
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




