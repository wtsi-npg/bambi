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
#include "htslib/hts.h"

/*
 * iterator structure
 */

typedef struct {
    samFile *f;
    bam_hdr_t *h;
    bam1_t *rec;
    bam1_t *nextRec;
} BAMit_t;

/*
 * Open a BAM file
 * arguments are: char *fname                filename to open
 *                char mode                  'r' or 'w'
 *                char *fmt                  format [bam,sam,cram] or NULL
 *                char compression level     [0..9] or NULL
 *                htsThreadPool *thread_pool thread pool to use, or NULL
 */
BAMit_t *BAMit_open(char *fname, char mode, char *fmt, char compression_level,
                    htsThreadPool *thread_pool);

/*
 * initialise with open file pointer and header
 */
BAMit_t *BAMit_init(samFile *f, bam_hdr_t *h);

/*
 * read next record and advance to next record
 */
bam1_t *BAMit_next(BAMit_t *bit);

/*
 * read next record *without* advancing to next record
 */

bam1_t *BAMit_peek(BAMit_t *bit);

/*
 * return true if there is a next record
 */

bool BAMit_hasnext(BAMit_t *bit);

/*
 * free data and close file
 */
void BAMit_free(void *bit);

#endif

