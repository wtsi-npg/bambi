/* filter.h

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

#ifndef __FILTERFILE_H__
#define __FILTERFILE_H__

#include <stdint.h>

typedef struct {
    int fhandle;
    char *errmsg;
    uint32_t version;
    uint32_t total_clusters;
    int current_cluster;
    int current_pf_cluster;
} filter_t;

filter_t *filter_open(char *fname);
int filter_next(filter_t *filter);
void filter_close(filter_t *filter);
void filter_seek(filter_t *filter, int cluster);

#endif

