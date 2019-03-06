/*  parse_bam.h - BAM parsing functions used by spatial_filter

    Copyright (C) 2017-2019 Genome Research Ltd.

    Author: Steven Leonard <srl@sanger.ac.uk>
            Jennifer Liddle <js10@sanger.ac.uk>

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

#ifndef __PARSE_BAM_H__
#define __PARSE_BAM_H__

#define BASE_ALIGN      (1<<0)
#define BASE_MISMATCH   (1<<1)
#define BASE_INSERTION  (1<<2)
#define BASE_DELETION   (1<<3)
#define BASE_SOFT_CLIP  (1<<4)
#define BASE_KNOWN_SNP  (1<<5)

#include <stdarg.h>

#include "hash_table.h"
#include "htslib/sam.h"
#include "bamit.h"

bam1_t *parse_bam_readinfo(BAMit_t *fp, 
                            int *bam_lane, 
                            int *bam_tile, 
                            int *bam_x, 
                            int *bam_y, 
                            int *bam_read, 
                         size_t *bam_offset);

int parse_bam_alignments(BAMit_t *fp, 
                            bam1_t *bam, 
                              char *read_seq, 
                               int *read_qual, 
                              char *read_ref, 
                               int *read_mismatch, 
                         const int read_buff_size,
                         HashTable *snp_hash);

void bam_header_add_pg(       char *id,
                              char *pn,
                              char *ds,
                              char *cl,
                      bam_hdr_t *bam_header);

char complement_base(char c);
void rev_comp_seq(char *seq);
char *reverse_seq(char *str);
int aux_type2size(uint8_t *s);
uint8_t *get_read(bam1_t *rec);
uint8_t *get_quality(bam1_t *rec);

#endif

