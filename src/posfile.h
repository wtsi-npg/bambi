/* posfile.h

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

#ifndef __POSFILE_H__
#define __POSFILE_H__

#include <stdint.h>

#define CLOCS_BLOCK_SIZE 25
#define CLOCS_IMAGE_WIDTH 2048
#define CLOCS_BLOCKS_PER_LINE ((CLOCS_IMAGE_WIDTH + CLOCS_BLOCK_SIZE - 1) / CLOCS_BLOCK_SIZE)

typedef enum { UNKNOWN_POS, POS, LOCS, CLOCS } POS_FILE_TYPE;

typedef struct {
    POS_FILE_TYPE file_type;
    int fhandle;
    char *errmsg;
    uint8_t version;
    uint32_t total_blocks;
    int current_block;
    uint8_t unread_clusters;
    int x,y;
} posfile_t;

posfile_t *posfile_open(char *fname);
int posfile_next(posfile_t *posfile);
void posfile_close(posfile_t *posfile);
void posfile_seek(posfile_t *posfile, int cluster);

static inline int posfile_get_x(posfile_t *posfile) { return posfile->x; }
static inline int posfile_get_y(posfile_t *posfile) { return posfile->y; }

#endif

