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

#include "bambi.h"
#include "filterfile.h"

filter_t *filter_open(char *fname)
{
    uint32_t empty;

    filter_t *filter = calloc(1, sizeof(filter_t));
    filter->total_clusters = 0;
    filter->current_cluster = 0;
    filter->buffer = NULL;
    filter->fhandle = open(fname,O_RDONLY);
    if (filter->fhandle == -1) {
        filter->errmsg = strdup(strerror(errno));
    } else {
        int n;
        filter->errmsg=NULL;
        n = read(filter->fhandle,(void *)&empty,4);
        n = read(filter->fhandle,(void *)&filter->version,4);
        n = read(filter->fhandle,(void *)&filter->total_clusters,4);
        if (n<0) {
            fprintf(stderr,"failed to read header from %s\n", fname);
            exit(1);
        }
    }
    return filter;
}

void filter_close(filter_t *filter)
{
    if (filter->fhandle!=-1) if (close(filter->fhandle)) { fprintf(stderr,"Can't close filter file\n"); exit(1); }
    free(filter->errmsg);
    free(filter->buffer);
    free(filter);
}

void filter_seek(filter_t *filter, int cluster)
{
    off_t pos = 12 + cluster;
    off_t n = lseek(filter->fhandle, pos, SEEK_SET);
    if (n != pos) {
        fprintf(stderr,"filter_seek(%d) failed: returned %d instead of %d\n", cluster, (int)n, (int)pos);
        exit(1);
    }
}

void filter_load(filter_t *filter, int64_t clusters)
{
    free(filter->buffer);
    filter->buffer = malloc(clusters);
    if (!filter->buffer) {
        fprintf(stderr, "filter_load(): Can't allocate %ld bytes for buffer\n", clusters);
        exit(1);
    }
    int64_t n = read(filter->fhandle, filter->buffer, clusters);
    if (n != clusters) {
        fprintf(stderr, "filter_load(): Expected %ld clusters, read %ld\n", clusters, n);
        exit(1);
    }
    filter->total_clusters = clusters;
}

char filter_get(filter_t *filter, int64_t n)
{
    return filter->buffer[n] & 0x01;
}

int filter_next(filter_t *filter)
{
    unsigned char next;

    if (read(filter->fhandle, (void *)&next, 1) != 1) {
        return -1;
    }

    filter->current_cluster++;
    return next & 0x01;
}

