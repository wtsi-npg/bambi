/*  rts.h - Region Table code used by spatial_filter

    Copyright (C) 2017 Genome Research Ltd.

    Author: Steven Leonard <srl@sanger.ac.uk>
            Jennifer Liddle <js10@sanger.ac.uk>

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

#ifndef RTS_H_INCLUDED
#define RTS_H_INCLUDED

#include <stdio.h>

#define N_READS 3
#define N_COMMENTS 100

#define COORD_SHIFT   1000
#define COORD_FACTOR  10

#define REGION_MAGIC                "RGFL"
#define REGION_SIZE                 200

#define REGION_STATE_COVERAGE   (1<<1)
#define REGION_STATE_MISMATCH   (1<<2)
#define REGION_STATE_INSERTION  (1<<3)
#define REGION_STATE_DELETION   (1<<4)
#define REGION_STATE_SOFT_CLIP  (1<<5)
#define REGION_STATE_BAD        (1<<6)

// The header of the filter file
typedef struct {
    char *region_magic;
    int coord_shift;
    int coord_factor;
    int ntiles;
    int ngood_tiles;
    int *tileArray;
    size_t *tileReadCountArray;
    int region_size;
    int nregions;
    int nregions_x;
    int nregions_y;
    int nreads;
    int readLength[N_READS];
    int totalReadLength;
    char *cmdLine;
    int ncomments;
    char *comments[N_COMMENTS];
    char *filterData;
    char *rgid;             // used as a key in the hash table
    uint64_t stats_nreads;       // number of records read
    uint64_t stats_nfiltered;    // number of records filtered out
} Header;

// An internal structure used to create the filter file
typedef struct {
        int align;
        int mismatch;
        int insertion;
        int deletion;
        int soft_clip;
        int known_snp;
        float quality;
	char state;
} RegionTable;

// Filter methods
void writeHeader(FILE *fp, Header *hdr);
void addHeaderComment(Header *hdr, char *comment);
void openFilters(va_t *filters, va_t *rgids);
//Header *readHeader(char *fname);
int setCurrentHdr(char *rgid);
char *getFilterData(int tile, int read, int cycle, int region);
Header *getHdr(char *rgid);
int x2region(int x, int region_size);
int xy2region(int x, int y);
int getHdrngood_tiles(void);
int getHdrReadLength(int read);
int getHdrnregions(void);
uint64_t getHdrStatsnreads(void);
void incHdrStatsnreads(void);
uint64_t getHdrStatsnfiltered(void);
void incHdrStatsnfiltered(void);
char *getHdrrgid(void);
va_t *HdrHash2Array(void);

#endif /* RTS_H_INCLUDED */
