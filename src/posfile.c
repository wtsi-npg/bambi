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

#include "bambi.h"
#include "posfile.h"

posfile_t *posfile_open(char *fname)
{
    posfile_t *posfile = calloc(1, sizeof(posfile_t));
    posfile->current_block = 0;
    posfile->errmsg = NULL;
    posfile->file_type = UNKNOWN_POS;
    posfile->fhandle = -1;
    posfile->x = NULL;
    posfile->y = NULL;

    // Should we handle compressed (.gz) files?

    char *base = basename(fname);
    char *ext = rindex(base,'.');
    if (ext) ext++;
    if (ext) {
        if (strcmp(ext,"clocs")==0) posfile->file_type = CLOCS;
        if (strcmp(ext,"locs")==0) posfile->file_type = LOCS;
        if (strcmp(ext,"txt")==0) posfile->file_type = POS;
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
        int n;
        n = read(posfile->fhandle,(void *)&posfile->version,1);
        n = read(posfile->fhandle,(void *)&posfile->total_blocks,4);
        n = read(posfile->fhandle,(void *)&posfile->unread_clusters,1);
        if (n<0) {
            fprintf(stderr,"failed to read header from %s\n", fname);
            exit(1);
        }
        posfile->current_block++;
    }

    if (posfile->file_type == LOCS) {
        int n;
        // first 8 bytes are unused
        n = read(posfile->fhandle,(void *)&posfile->total_blocks,4);
        n = read(posfile->fhandle,(void *)&posfile->total_blocks,4);
        n = read(posfile->fhandle,(void *)&posfile->total_blocks,4);
        if (n<0) {
            fprintf(stderr,"failed to read header from %s\n", fname);
            exit(1);
        }
    }
    return posfile;
}

/*
 * seek to a given cluster number
 */
void posfile_seek(posfile_t *posfile, int cluster)
{
    off_t pos = 12 + cluster * 8;
    if (posfile->file_type != LOCS) {
        fprintf(stderr,"Can only handle NextSeq pos files of type LOC\n");
        exit(1);
    }

    off_t r = lseek(posfile->fhandle, pos, SEEK_SET);
    if (r != pos) {
        fprintf(stderr,"Trying to seek to %d (cluster %d) but returned %d\n", (int)pos, cluster, (int)r);
        perror("posfile_seek() failed");
    }
}

void posfile_close(posfile_t *posfile)
{
    free(posfile->errmsg);
    free(posfile->x); free(posfile->y);
    if (posfile->fhandle >= 0) {
        if(close(posfile->fhandle)) { fprintf(stderr,"Can't close posfile"); exit(1); }
    }
    free(posfile);
}

static void locs_load(posfile_t *posfile, filter_t *filter)
{
    float dx, dy;
    int i,j,f;
    int bufsize = posfile->total_blocks * 4 * 2;
    char *buffer = malloc(bufsize);

    free(posfile->x); posfile->x = malloc(posfile->total_blocks * sizeof(int));
    free(posfile->y); posfile->y = malloc(posfile->total_blocks * sizeof(int));
    if (!buffer || !posfile->x || !posfile->y) {
        fprintf(stderr,"locs_load(): failed to malloc buffer for %d blocks\n", posfile->total_blocks);
        exit(1);
    }

    int n = read(posfile->fhandle, buffer, bufsize);
    if (n != bufsize) {
        fprintf(stderr,"locs_load(): expected %d, read %d\n", bufsize, n);
        exit(1);
    }

    j=0;
    for (i=0, f=0; i < bufsize; i+=8, f++) {
        dx = *(float *)(buffer+i);
        dy = *(float *)(buffer+i+4);
        if (filter && f >= filter->total_clusters) break;
        if (!filter || (filter->buffer[f] & 0x01)) {
            posfile->x[j] = 10 * dx + 1000.5;
            posfile->y[j] = 10 * dy + 1000.5;
            j++;
        }
    }

    posfile->size = j;
    free(buffer);
}

void clocs_load(posfile_t *posfile, int bufsize, filter_t *filter)
{
    unsigned char dx, dy;
    int j = 0;
    int f = 0;
    if (bufsize == 0) bufsize = 10000;
    free(posfile->x); posfile->x = malloc(bufsize * sizeof(int));
    free(posfile->y); posfile->y = malloc(bufsize * sizeof(int));
    if (!posfile->x || !posfile->y) {
        fprintf(stderr,"clocs_load(): failed to malloc buffer for %d blocks\n", posfile->total_blocks);
        exit(1);
    }

    for (;;) {
        while (posfile->unread_clusters == 0 && (posfile->current_block < posfile->total_blocks)) {
            if (read(posfile->fhandle, (void *)&posfile->unread_clusters, 1) != 1) break;
            posfile->current_block++;
        }

        if (posfile->unread_clusters == 0) break;
        posfile->unread_clusters--;

        int n;
        n = read(posfile->fhandle, (void *)&dx, 1);    
        n = read(posfile->fhandle, (void *)&dy, 1);
        if (n<0) {
            fprintf(stderr,"something has gone wrong in clocs_load()\n");
            exit(1);
        }

        if (j >= bufsize) {
            bufsize *= 2;
            posfile->x = realloc(posfile->x,bufsize * sizeof(int));
            posfile->y = realloc(posfile->y,bufsize * sizeof(int));
            if (!posfile->x || !posfile->y) {
                fprintf(stderr,"clocs_load(): failed to realloc %d entries\n", bufsize);
                exit(1);
            }
        }
        if (!filter || (filter->buffer[f] & 0x01)) {
            posfile->x[j] = 10 * CLOCS_BLOCK_SIZE * ((posfile->current_block - 1) % CLOCS_BLOCKS_PER_LINE) + dx + 1000;
            posfile->y[j] = 10 * CLOCS_BLOCK_SIZE * ((posfile->current_block - 1) / CLOCS_BLOCKS_PER_LINE) + dy + 1000;
            j++;
        }
        f++;
    }
    posfile->size = j;
}

/*
int posfile_next(posfile_t *posfile)
{
    if (posfile->file_type == CLOCS) return clocs_next(posfile);
    if (posfile->file_type == LOCS) return locs_next(posfile);
    return -1;
}
*/

void posfile_load(posfile_t *posfile, int bufsize, filter_t *filter)
{
    if (posfile->file_type == CLOCS) return clocs_load(posfile, bufsize, filter);
    if (posfile->file_type == LOCS) return locs_load(posfile, filter);
}

