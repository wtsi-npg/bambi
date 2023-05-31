/* posfile.c

   Functions to read and parse an Illumina pos, locs, or clocs file.

    Copyright (C) 2016-2023 Genome Research Ltd.

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
    if (!posfile) {
        fprintf(stderr, "Out of memory");
        exit(-1);
    }
    posfile->current_block = 0;
    posfile->file_name = strdup(fname);
    posfile->errmsg = NULL;
    posfile->file_type = UNKNOWN_POS;
    posfile->fhandle = NULL;
    posfile->x = NULL;
    posfile->y = NULL;

    if (!posfile->file_name) {
        fprintf(stderr, "Out of memory");
        exit(-1);
    }

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
    posfile->fhandle = fopen(fname, "rb");

    if (posfile->fhandle == NULL) {
        posfile->errmsg = strdup(strerror(errno));
        return posfile;
    }

    if (posfile->file_type == CLOCS) {
        int n;
        n = fread(&posfile->version,1, 1, posfile->fhandle);
        if (n == 1) n = fread(&posfile->total_blocks, 4, 1, posfile->fhandle);
        if (n == 1) n = fread(&posfile->unread_clusters,1, 1, posfile->fhandle);
        if (n != 1) {
            fprintf(stderr,"failed to read header from %s\n", fname);
            exit(1);
        }
        posfile->current_block++;
    }

    if (posfile->file_type == LOCS) {
        int n;
        uint32_t x[3];
        // first 8 bytes are unused
        n = fread(x, 4, 3, posfile->fhandle);
        if (n != 3) {
            fprintf(stderr,"failed to read header from %s\n", fname);
            exit(1);
        }
        posfile->total_blocks = x[2];
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

    int r = fseeko(posfile->fhandle, pos, SEEK_SET);
    if (r < 0) {
        fprintf(stderr,"Trying to seek on %s to %ld (cluster %d) but returned %d\n", posfile->file_name, (long) pos, cluster, r);
        perror("posfile_seek() failed");
        exit(1);
    }
}

void posfile_close(posfile_t *posfile)
{
    free(posfile->errmsg);
    free(posfile->x); free(posfile->y);
    if (posfile->fhandle) {
        if(fclose(posfile->fhandle)) {
            fprintf(stderr,"Can't close posfile %s : %s", posfile->file_name, strerror(errno));
            exit(1);
        }
    }
    free(posfile->file_name);
    free(posfile);
}

static void locs_load(posfile_t *posfile, filter_t *filter)
{
    float dx, dy;
    int i,j,f;
    size_t bufsize = posfile->total_blocks * 4 * 2;
    char *buffer = malloc(bufsize);

    free(posfile->x); posfile->x = malloc(posfile->total_blocks * sizeof(int));
    free(posfile->y); posfile->y = malloc(posfile->total_blocks * sizeof(int));
    if (!buffer || !posfile->x || !posfile->y) {
        fprintf(stderr,"locs_load(): failed to malloc buffer for %d blocks\n", posfile->total_blocks);
        exit(1);
    }

    size_t n = fread(buffer, 1, bufsize, posfile->fhandle);
    if (n != bufsize) {
        fprintf(stderr,"locs_load(%s): expected %zd, read %zd\n", posfile->file_name, bufsize, n);
        exit(1);
    }

    j=0;
    for (i=0, f=0; i < bufsize; i+=8, f++) {
        dx = *(float *)(buffer+i);
        dy = *(float *)(buffer+i+4);
        if (filter && f >= filter->total_clusters) break;
        if (!filter || (filter->buffer[f] & 0x01)) {
            // leaving dx,dy as floats would lose us precision, so force to doubles
            posfile->x[j] = 10 * (double)dx + 1000.5;
            posfile->y[j] = 10 * (double)dy + 1000.5;
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
            if (fread(&posfile->unread_clusters, 1, 1, posfile->fhandle) != 1) break;
            posfile->current_block++;
        }

        if (posfile->unread_clusters == 0) break;
        posfile->unread_clusters--;

        size_t n;
        n = fread(&dx, 1, 1, posfile->fhandle);
        if (n == 1) n = fread(&dy, 1, 1, posfile->fhandle);
        if (n != 1) {
            if (ferror(posfile->fhandle)) {
                fprintf(stderr,"clocs_load(%s): %s\n", posfile->file_name, strerror(errno));
                exit(1);
            } else {
                fprintf(stderr,"clocs_load(%s): Warning: reached end of file with %u clusters and %u blocks unread\n",
                        posfile->file_name, (int) posfile->unread_clusters,
                        posfile->total_blocks - posfile->current_block);
                break;
            }
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

