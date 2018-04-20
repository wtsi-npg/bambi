/*  rts.c - Region Table code used by spatial_filter

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

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <search.h>

#include "bambi.h"
#include "rts.h"
#include "hash_table.h"
#include "array.h"

static Header *_hdr = NULL;
static HashTable *_filter_hash = NULL;

// The Perl 'chomp()' function is too useful not to recreate it here...
static void chomp(char *line)
{
    int n = strlen(line) - 1;
    if (line[n] == '\n') line[n] = 0;
    return;
}

inline int getHdrReadLength(int read) { return _hdr->readLength[read]; }
inline int getHdrnregions(void) { return _hdr->nregions; }
//
// Filter methods
//

// Allocate memory for the filter file, then slurp the whole thing
static void readFilterData(FILE *fp, Header *hdr)
{
    int i, n=0;
    for (i=0; i < hdr->nreads; i++) n += hdr->readLength[i];
    if ((hdr->ntiles * n * hdr->nregions)) {
        hdr->filterData = malloc(hdr->ntiles * n * hdr->nregions);
        if (fread(hdr->filterData, hdr->ntiles * n * hdr->nregions, 1, fp) != 1) {
            die("Error reading filter file\n");
        }
    }
}


Header *readHeader(char *fname)
{
    int len = 1024;
    char line[1024];
    int i;
    char *p;
    Header *hdr;
    FILE *fp;

    hdr = malloc(sizeof(Header));
    if (!hdr) die("readHeader(): Can't allocate memory for %s\n", fname);

    fp = fopen(fname, "rb");
    if (!fp) die("readHeader(): Can't open %s\n", fname);

    p=fgets(line, len, fp); chomp(line); hdr->region_magic = strdup(line);
    p=fgets(line, len, fp); hdr->coord_shift = atoi(line);
    p=fgets(line, len, fp); hdr->coord_factor = atoi(line);
    p=fgets(line, len, fp); hdr->region_size = atoi(line);
    p=fgets(line, len, fp); hdr->ntiles = atoi(line);
    hdr->tileArray = NULL;
    hdr->tileReadCountArray = NULL;
    if (hdr->ntiles > 0) {
        hdr->tileArray = malloc(hdr->ntiles * sizeof(int));
        hdr->tileReadCountArray = malloc(hdr->ntiles * sizeof(size_t));
        for (i=0; i < hdr->ntiles; i++) {
            int n;
            p=fgets(line, len, fp);
            n = sscanf(line, "%d\t%lu\n", &hdr->tileArray[i], &hdr->tileReadCountArray[i]);
            switch (n) {
                case 1:	hdr->tileReadCountArray[i] = 0; break;
                case 2:	break;
                default: fprintf(stderr,"ERROR: Invalid filter file");
                         exit(1);
                         break;
            }
        }
    }
    p=fgets(line, len, fp); hdr->nregions = atoi(line);
    p=fgets(line, len, fp); hdr->nregions_x = atoi(line);
    p=fgets(line, len, fp); hdr->nregions_y = atoi(line);
    p=fgets(line, len, fp); hdr->nreads = atoi(line);
    hdr->totalReadLength = 0;
    for (i=0; i < hdr->nreads; i++) {
        p=fgets(line, len, fp); hdr->readLength[i] = atoi(line);
        hdr->totalReadLength += hdr->readLength[i];
    }
    p=fgets(line, len, fp); chomp(line); hdr->cmdLine = strdup(line);
    p=fgets(line, len, fp); hdr->ncomments = atoi(line);
    for (i=0; i < hdr->ncomments; i++) {
        p=fgets(line, len, fp); chomp(line); hdr->comments[i] = strdup(line);
    }

    readFilterData(fp, hdr);

    fclose(fp);
    return hdr;
}

void writeHeader(FILE *fp, Header *hdr)
{
    int i;
    fprintf(fp, "%s\n", hdr->region_magic);
    fprintf(fp, "%d\n", hdr->coord_shift);
    fprintf(fp, "%d\n", hdr->coord_factor);
    fprintf(fp, "%d\n", hdr->region_size);
    fprintf(fp, "%d\n", hdr->ntiles);
    for (i=0; i < hdr->ntiles; i++) {
      fprintf(fp, "%d\t%lu\n", hdr->tileArray[i], hdr->tileReadCountArray[i]);
	}
    fprintf(fp, "%d\n", hdr->nregions);
    fprintf(fp, "%d\n", hdr->nregions_x);
    fprintf(fp, "%d\n", hdr->nregions_y);
    fprintf(fp, "%d\n", hdr->nreads);
    for (i=0; i<hdr->nreads; i++)
        fprintf(fp, "%d\n", hdr->readLength[i]);
    fprintf(fp, "%s\n", hdr->cmdLine);
    fprintf(fp, "%d\n", hdr->ncomments);
    for (i=0; i<hdr->ncomments; i++)
        fprintf(fp, "%s\n", hdr->comments[i]);
}

void addHeaderComment(Header *hdr, char *comment)
{
    hdr->comments[hdr->ncomments++] = strdup(comment);
}

static int keyComp(const void *k1, const void *k2)
{
	int key = *(int *)k1;
	int mem = *(int *)k2;
	return key != mem;
}

char* getFilterData(int tile, int read, int cycle, int region)
{
    int itile, i, previousReadLength = 0, offset;
    size_t nelem = _hdr->ntiles;
    void *pitile = lfind(&tile, _hdr->tileArray, &nelem, sizeof(int), &keyComp);
    if (!pitile) return NULL;  	// if tile not found in filter
    itile = ((int*)pitile - _hdr->tileArray);
    for (i=0; i < read; i++) previousReadLength += _hdr->readLength[i];
	offset = itile * _hdr->totalReadLength * _hdr->nregions + (previousReadLength + cycle) * _hdr->nregions + region;
	return _hdr->filterData + offset;
}

// which region is x in?
int inline x2region(int x, int region_size)
{
    int coord_shift = (_hdr ? _hdr->coord_shift : COORD_SHIFT);
    int coord_factor = (_hdr ? _hdr->coord_factor : COORD_FACTOR);
    float x_coord = (float)(x - coord_shift) / (float)coord_factor;
    return (int)(x_coord / region_size);
}

// which region is (x,y) in?
inline int xy2region(int x, int y)
{
    return x2region(x, _hdr->region_size) * _hdr->nregions_y + x2region(y, _hdr->region_size);
}

Header *getHdr(char *rgid)
{
    if (!rgid) rgid = "null";
    HashItem *hi = HashTableSearch(_filter_hash, rgid, 0);
    if (hi) return hi->data.p;
    return NULL;
}

int setCurrentHdr(char *rgid)
{
    _hdr = getHdr(rgid);
    return _hdr != NULL;
}

void openFilters(va_t *fnames, va_t *rgids)
{
    _filter_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS);

    for (int n = 0; n < fnames->end; n++) {
        HashData hd;
        char *fname = fnames->entries[n];
        char *rgid = rgids ? rgids->entries[n] : "null";
        Header *hdr = readHeader(fname);
        hd.p = hdr;
        HashTableAdd(_filter_hash, rgid, 0, hd, NULL);
    }
}


