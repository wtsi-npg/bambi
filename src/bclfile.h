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

typedef enum { BCL, SCL } BCL_FILE_TYPE;

typedef struct {
    BCL_FILE_TYPE file_type;
    int fhandle;
    gzFile gzhandle;
    char *errmsg;
    uint32_t total_clusters;
    int current_cluster;
    int current_base;
    char base;
    int quality;
char *filename;
} bclfile_t;

bclfile_t *bclfile_open(char *fname);
int bclfile_next(bclfile_t *bclfile);
void bclfile_close(bclfile_t *bclfile);
void bclfile_seek(bclfile_t *bclfile, int cluster);

#endif

