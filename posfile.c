/* posfile.c

   Functions to read and parse an Illumina pos, locs, or clocs file.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <libgen.h>

#include "posfile.h"

posfile_t *posfile_open(char *fname)
{
    posfile_t *posfile = calloc(1, sizeof(posfile_t));
    posfile->current_block = 0;
    posfile->errmsg = NULL;
    posfile->file_type = UNKNOWN_POS;

    // Should we handle compressed (.gz) files?

    char *base = basename(fname);
    char *ext = rindex(base,'.');
    if (ext) ext++;
    if (ext) {
        if (strcmp(ext,"clocs")==0) posfile->file_type = CLOCS;
        if (strcmp(ext,"locs")==0) posfile->file_type = LOCS;
    }
    if (posfile->file_type == UNKNOWN_POS) {
        posfile->errmsg = strdup("posfile_open(): Unknown file type\n");
        return posfile;
    }
    posfile->fhandle = open(fname,O_RDONLY);

    if (posfile->fhandle == -1) {
        posfile->errmsg = strdup(strerror(errno));
        return posfile;
    }

    if (posfile->file_type == CLOCS) {
        read(posfile->fhandle,(void *)&posfile->version,1);
        read(posfile->fhandle,(void *)&posfile->total_blocks,4);
        read(posfile->fhandle,(void *)&posfile->unread_clusters,1);
        posfile->current_block++;
    }

    if (posfile->file_type == LOCS) {
        // first 8 bytes are unused
        read(posfile->fhandle,(void *)&posfile->total_blocks,4);
        read(posfile->fhandle,(void *)&posfile->total_blocks,4);
        read(posfile->fhandle,(void *)&posfile->total_blocks,4);
    }

    return posfile;
}

void posfile_close(posfile_t *posfile)
{
    free(posfile->errmsg);
    if (posfile->fhandle >= 0) close(posfile->fhandle);
    free(posfile);
}

static int locs_next(posfile_t *posfile)
{
    float dx, dy;

    if (posfile->current_block >= posfile->total_blocks) return -1;
    posfile->current_block++;

    read(posfile->fhandle, (void *)&dx, 4);
    read(posfile->fhandle, (void *)&dy, 4);

    posfile->x = 10 * dx + 1000.5;
    posfile->y = 10 * dy + 1000.5;

    return 0;
}

static int clocs_next(posfile_t *posfile)
{
    unsigned char dx, dy;

    while (posfile->unread_clusters == 0 && (posfile->current_block < posfile->total_blocks)) {
        if (read(posfile->fhandle, (void *)&posfile->unread_clusters, 1) != 1) return -1;
        posfile->current_block++;
    }

    if (posfile->unread_clusters == 0) return -1;
    posfile->unread_clusters--;

    read(posfile->fhandle, (void *)&dx, 1);    
    read(posfile->fhandle, (void *)&dy, 1);

    posfile->x = 10 * CLOCS_BLOCK_SIZE * ((posfile->current_block - 1) % CLOCS_BLOCKS_PER_LINE) + dx + 1000;
    posfile->y = 10 * CLOCS_BLOCK_SIZE * ((posfile->current_block - 1) / CLOCS_BLOCKS_PER_LINE) + dy + 1000;
    return 0;
}

int posfile_next(posfile_t *posfile)
{
    if (posfile->file_type == CLOCS) return clocs_next(posfile);
    if (posfile->file_type == LOCS) return locs_next(posfile);
    return -1;
}

