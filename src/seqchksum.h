/* seqchksum.h

    Copyright (C) 2018 Genome Research Ltd.

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

#ifndef __SEQCHKSUM_H__
#define __SEQCHKSUM_H__

#include "hash_table.h"

typedef enum {
    HASH_UNKNOWN,
    HASH_CRC32,
    HASH_CRC32PROD
} HASH_TYPE;

#define DEFAULT_HASH_TYPE   HASH_CRC32PROD

/* 
 * digest line struct. Element 0 is 'all', element 1 is 'pass'
 */
typedef struct {
    uint32_t count[2];
    uint32_t chksum[4][2];  // first index is column, second is 'all','pass'
} digest_line_t;

typedef struct {
    digest_line_t all;
    HashTable *rgHash;
} chksum_results_t;

/*
 * Initialise results structure
 */
chksum_results_t *chksum_init_results(HASH_TYPE hash);

/*
 * free results structure
 */
void chksum_free_results(chksum_results_t *results);

/*
 * process one BAM record, and store accumulated results in 'results'
 */
int seqchksum_processRecord(bam1_t *rec, HASH_TYPE hash, chksum_results_t *results);

#endif

