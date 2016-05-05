/* filterfile.c

   Functions to read and parse an Illumina filter file.

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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#include "filterfile.h"

filter_t *filter_open(char *fname)
{
    uint32_t empty;

    filter_t *filter = calloc(1, sizeof(filter_t));
    filter->total_clusters = 0;
    filter->current_cluster = 0;
    filter->current_pf_cluster = 0;
    filter->fhandle = open(fname,O_RDONLY);
    if (filter->fhandle == -1) {
        filter->errmsg = strdup(strerror(errno));
    } else {
        read(filter->fhandle,(void *)&empty,4);
        read(filter->fhandle,(void *)&filter->version,4);
        read(filter->fhandle,(void *)&filter->total_clusters,4);
    }
    return filter;
}

void filter_close(filter_t *filter)
{
    close(filter->fhandle);
}

int filter_next(filter_t *filter)
{
    unsigned char next;

    if (read(filter->fhandle, (void *)&next, 1) != 1) {
        return -1;
    }

    filter->current_cluster++;
    next = next & 0x01;
    if (next == 1) filter->current_pf_cluster++;
    return next;
}

