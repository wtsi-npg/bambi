/*  bamit.h -- BAM Iterator - for reading BAM files one record at a time.

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

#ifndef __BAMIT_H__
#define __BAMIT_H__

#include <stdbool.h>
#include <stdlib.h>
#include "htslib/sam.h"

/*
 * iterator structure
 */

typedef struct {
    samFile *f;
    bam_hdr_t *h;
    bam1_t *rec;
    bam1_t *nextRec;
} BAMit_t;

BAMit_t *BAMit_init(samFile *f, bam_hdr_t *h);
void BAMit_free(void *bit);
bam1_t *BAMit_next(BAMit_t *bit);
bam1_t *BAMit_peek(BAMit_t *bit);
bool BAMit_hasnext(BAMit_t *bit);

#endif

