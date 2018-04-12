/* filterfile.c

   Functions to read and parse an Illumina filter file.

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
#include <assert.h>

#include "bambi.h"
#include "filterfile.h"

filter_t *filter_open(char *fname)
{
    filter_t *filter = calloc(1, sizeof(filter_t));
    if (!filter) {
        fprintf(stderr, "Out of memory");
        exit(1);
    }
    filter->total_clusters = 0;
    filter->current_cluster = 0;
    filter->buffer = NULL;
    filter->buffer_size = 0;
    filter->fhandle = fopen(fname, "rb");
    if (filter->fhandle == NULL) {
        filter->errmsg = strdup(strerror(errno));
    } else {
        int n;
        uint32_t x[3];
        filter->errmsg=NULL;
        n = fread(x, 4, 3, filter->fhandle);
        if (n != 3) {
            fprintf(stderr,"failed to read header from %s\n", fname);
            exit(1);
        }
        // x[0] is ignored
        filter->version = x[1];
        filter->total_clusters = x[2];
    }
    return filter;
}

void filter_close(filter_t *filter)
{
    if (filter->fhandle!=NULL) {
        if (fclose(filter->fhandle)) {
            fprintf(stderr, "Can't close filter file\n");
            exit(1);
        }
    }
    free(filter->errmsg);
    free(filter->buffer);
    free(filter);
}

void filter_seek(filter_t *filter, int cluster)
{
    off_t pos = 12 + cluster;
    int n = fseeko(filter->fhandle, pos, SEEK_SET);
    if (n < 0) {
        fprintf(stderr, "filter_seek(%d) failed: %s\n", cluster, strerror(errno));
        exit(1);
    }
}

void filter_load(filter_t *filter, size_t clusters)
{
    free(filter->buffer);
    filter->buffer = malloc(clusters);
    if (!filter->buffer) {
        fprintf(stderr, "filter_load(): Can't allocate %zd bytes for buffer\n", clusters);
        exit(1);
    }
    size_t n = fread(filter->buffer, 1, clusters, filter->fhandle);
    if (n != clusters) {
        fprintf(stderr, "filter_load(): Expected %ld clusters, read %zd\n", clusters, n);
        exit(1);
    }
    filter->total_clusters = clusters;
    filter->buffer_size = clusters;
}

char filter_get(filter_t *filter, size_t n)
{
    assert(n < filter->buffer_size);
    return filter->buffer[n] & 0x01;
}

int filter_next(filter_t *filter)
{
    unsigned char next;

    if (fread(&next, 1, 1, filter->fhandle) != 1) {
        return -1;
    }

    filter->current_cluster++;
    return next & 0x01;
}

