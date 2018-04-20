/* bclfile.h

    Copyright (C) 2016 Genome Research Ltd.

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

#ifndef __BCLFILE_H__
#define __BCLFILE_H__

#include <stdint.h>
#include <zlib.h>
#include "array.h"
#include "filterfile.h"
#include "bambi.h"

typedef struct {
    uint32_t  tilenum;
    uint32_t  nclusters;
    uint32_t  uncompressed_blocksize;
    uint32_t  compressed_blocksize;
} tilerec_t;
    
typedef struct {
    MACHINE_TYPE machine_type;
    FILE *fhandle;
    gzFile gzhandle;
    char *errmsg;
    int is_cached;
    uint32_t total_clusters;
    int current_base;

    int bases_size;
    int base_ptr;
    char *bases;
    char *quals;
    char base;
    int quality;
    char *filename;

    char current_byte;
    int block_index;
    // CBCL specific fields
    uint16_t version;
    uint32_t header_size;
    unsigned char bits_per_base;
    unsigned char bits_per_qual;
    uint32_t nbins;
    int qbin[4];
    uint32_t ntiles;
    tilerec_t *current_tile;
    va_t *tiles;
    char *current_block;
    char *current_block_ptr;
    uint32_t current_block_size;
    char pfFlag;
    int surface;
    int fails;
} bclfile_t;

int bcl_tile2surface(int tile);
bclfile_t *bclfile_open(char *fname, MACHINE_TYPE mt);
void bclfile_close(bclfile_t *bclfile);
int bclfile_load_tile(bclfile_t *bclfile, int tile, filter_t *filter);
char bclfile_base(bclfile_t *bcl, int cluster);
int bclfile_quality(bclfile_t *bcl, int cluster);
#endif

