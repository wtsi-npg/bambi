/* posfile.h

    Copyright (C) 2016 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef __POSFILE_H__
#define __POSFILE_H__

#include <stdint.h>

#define CLOCS_BLOCK_SIZE 25
#define CLOCS_IMAGE_WIDTH 2048
#define CLOCS_BLOCKS_PER_LINE ((CLOCS_IMAGE_WIDTH + CLOCS_BLOCK_SIZE - 1) / CLOCS_BLOCK_SIZE)

typedef enum { POS, LOCS, CLOCS } FILE_TYPE;

typedef struct {
    FILE_TYPE file_type;
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

static inline int posfile_get_x(posfile_t *posfile) { return posfile->x; }
static inline int posfile_get_y(posfile_t *posfile) { return posfile->y; }

#endif

