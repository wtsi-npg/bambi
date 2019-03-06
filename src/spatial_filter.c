/*  spatial_filter.c - This code looks for spatial features given an aligned bam file

    Copyright (C) 2018 Genome Research Ltd.

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

#include "bambi.h"
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <search.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <fcntl.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdarg.h>
#include <getopt.h>
#include <gd.h>
#include <gdfonts.h>
#include <limits.h>

#include "htslib/sam.h"
#include "htslib/hfile.h"
#include "bamit.h"
#include "hash_table.h"
#include "parse_bam.h"
#include "array.h"
#include "parse.h"

#define SF_MAX_LANES 17

#define N_READS 3

#define COORD_SHIFT   1000
#define COORD_FACTOR  10

#define REGION_MAGIC                "RGF3"
#define REGION_MAGIC_LEN            5
#define SF_CMDLINE_LEN              1024

#define REGION_SIZE                 200

#define REGION_STATE_COVERAGE   (1<<1)
#define REGION_STATE_MISMATCH   (1<<2)
#define REGION_STATE_INSERTION  (1<<3)
#define REGION_STATE_DELETION   (1<<4)
#define REGION_STATE_SOFT_CLIP  (1<<5)
#define REGION_STATE_BAD        (1<<6)

#define REGION_MISMATCH_THRESHOLD   0.016  // threshold for setting region mismatch state
#define REGION_INSERTION_THRESHOLD  0.016  // threshold for setting region insertion state
#define REGION_DELETION_THRESHOLD   0.016  // threshold for setting region deletion state

#define TILE_REGION_THRESHOLD  0.75  // threshold for setting region state at tile level

#define MIN_TILE_READ_COUNT  1000 // min number of aligned reads on a tile

#define REGION_STATE_MASK  (REGION_STATE_INSERTION | REGION_STATE_DELETION)  // region mask used to filter reads

// key for the region hash table
typedef struct {
    uint32_t x;
    uint32_t y;
} region_key_t;

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
} RegionTableEntry_t, ***RegionTable_t;

// The header of the filter file
typedef struct {
    char region_magic[REGION_MAGIC_LEN];
    char cmdLine[SF_CMDLINE_LEN];
} filter_header_t;

// Header for each Lane of the filter file
typedef struct {
    int lane;
    int coord_shift;
    int coord_factor;
    size_t ntiles;
    int *tileArray;
    size_t *tileReadCountArray;
    HashTable *region_hash;
    int nregions;
    int *regions;
    int region_size;
    int nregions_x;
    int nregions_y;
    uint64_t nreads;
    int readLength[N_READS];
    int totalReadLength;
    uint32_t filterDataSize;
    char *filterData;
    uint64_t stats_nreads;       // number of records read
    uint64_t stats_nfiltered;    // number of records filtered out
    int ngood_tiles;
    RegionTable_t rts;
} Header;

static Header *LaneArray[SF_MAX_LANES];
static filter_header_t Fheader;

enum images { IMAGE_COVERAGE,
              IMAGE_DELETION,
              IMAGE_INSERTION,
              IMAGE_MISMATCH,
              IMAGE_QUALITY,
              N_IMAGES };

char *image_names[] = {"cov",
                       "del",
                       "ins",
                       "mma",
                       "qua"};

#define IMAGE_COLUMN_GAP 3
#define IMAGE_LABEL_HEIGHT 25
#define NUM_IMAGES_IN_REPORT_ROW 18
    
enum colours { COLOUR_LEVEL_0,
               COLOUR_LEVEL_1,
               COLOUR_LEVEL_2,
               COLOUR_LEVEL_3,
               COLOUR_LEVEL_4,
               COLOUR_LEVEL_5,
               COLOUR_LEVEL_6,
               COLOUR_LEVEL_7,
               COLOUR_LEVEL_8,
               COLOUR_LEVEL_9,
               COLOUR_LEVEL_10,
               COLOUR_LEVEL_11,
               COLOUR_TEXT,
               COLOUR_QC_FAIL,
               COLOUR_ZERO_QUAL,
               COLOUR_LOW_QUAL,
               COLOUR_MEDIUM_QUAL,
               COLOUR_HIGH_QUAL,
               N_COLOURS };

static int *colour_table = NULL;

typedef struct {
	va_t *filters;
	char *snp_file;
	char *in_bam_file;
	HashTable *snp_hash;
	char *working_dir;
	char *output;
	char *apply_stats_out;
	int calculate;
    bool dumpFilter;
	char *tileviz;
	int apply;
	int qcfail;
	int verbose;
	int region_min_count;
    int region_size;
	float region_mismatch_threshold;
	float region_insertion_threshold;
	float region_deletion_threshold;
	char compression_level;
    char *argv_list;
    char *input_fmt;
    char *output_fmt;
} opts_t;

#define min(a, b) ( (a<=b) ? a : b )
#define max(a, b) ( (a>=b) ? a : b )

int int_cmp(const void *i1, const void *i2) {
    return *(int *)i1 != *(int *)i2;
}

int int_sort(const void *i1, const void *i2) {
    return *(int *)i1 > *(int *)i2;
}

/*
 * Create and Initialise a header structure with default values
 */
static Header *sf_initHdr(void)
{
    Header *hdr = smalloc(sizeof(Header));
    hdr->lane = 0;
    hdr->coord_shift = COORD_SHIFT;
    hdr->coord_factor = COORD_FACTOR;
    hdr->ntiles = 0;
    hdr->ngood_tiles = 0;
    hdr->tileArray = NULL;
    hdr->tileReadCountArray = NULL;
    hdr->region_size = REGION_SIZE;
    hdr->region_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE|HASH_FUNC_JENKINS);
    hdr->regions = NULL;
    hdr->nregions = 0;
    hdr->nregions_x = 0;
    hdr->nregions_y = 0;
    hdr->nreads = 0;
    for (int n=0; n < N_READS; n++) hdr->readLength[n] = 0;
    hdr->totalReadLength = 0;
    hdr->filterDataSize = 0;
    hdr->filterData = NULL;
    hdr->stats_nreads = 0;
    hdr->stats_nfiltered = 0;

    return hdr;
}

// which region is x in?
static uint32_t x2region(Header *hdr, int x)
{
    float x_coord = (float)(x - hdr->coord_shift) / (float)hdr->coord_factor;
    return (int)(x_coord / hdr->region_size);
}

// which region is (x,y) in?
static uint32_t xy2region(Header *hdr, int x, int y)
{
    return x2region(hdr, x) * hdr->nregions_y + x2region(hdr, y);
}

// Convert tile to an array index
static int tile2index(Header *hdr, int tile)
{
    void *pitile = lfind(&tile, hdr->tileArray, &hdr->ntiles, sizeof(int), &int_cmp);
    if (!pitile) return -1;
    return ((int*)pitile - hdr->tileArray);
}

// Return filter data for given lane, tile, read, cycle, region
static char *getFilterData(Header *hdr, int itile, int read, int cycle, int region)
{
    int previousReadLength = 0, offset;
    for (int n=0; n < read; n++) previousReadLength += hdr->readLength[n];
    offset = itile * hdr->totalReadLength * hdr->nregions + (previousReadLength + cycle) * hdr->nregions + region;
    return hdr->filterData + offset;
}

/*
 * Read the filter header from a filter file
 */
static void readFheader(hFILE *fp)
{
    if (hread(fp, &Fheader, sizeof(Fheader)) < 0) die("readFheader() failed\n");
    if (strncmp(Fheader.region_magic,REGION_MAGIC, strlen(REGION_MAGIC)-1)) die("Not a valid filter file\n");
}

/*
 * Read a Lane header from a filter file
 */
static Header *readHeader(hFILE *fp)
{
    Header *hdr;

    hdr = smalloc(sizeof(Header));
    hdr->ngood_tiles = 0;

    ssize_t r = hread(fp, &hdr->lane, sizeof(hdr->lane));
    if (r==0) { free(hdr); return NULL; }   // end of file, no more lanes left
    if (r <0) goto fail;

    if (hread(fp, &hdr->coord_shift, sizeof(hdr->coord_shift)) < 0) goto fail;
    if (hread(fp, &hdr->coord_factor, sizeof(hdr->coord_factor)) < 0) goto fail;
    if (hread(fp, &hdr->ntiles, sizeof(hdr->ntiles)) < 0) goto fail;

    if (hdr->ntiles > 0) {
        hdr->tileArray = smalloc(hdr->ntiles * sizeof(int));
        hdr->tileReadCountArray = malloc(hdr->ntiles * sizeof(size_t));
        for (int i=0; i < hdr->ntiles; i++) {
            if (hread(fp, &hdr->tileArray[i], sizeof(hdr->tileArray[i])) < 0) die("readHeader() failed\n");
            if (hread(fp, &hdr->tileReadCountArray[i], sizeof(hdr->tileReadCountArray[i])) < 0) die("readHeader() failed\n");
        }
    }
    if (hread(fp, &hdr->nregions, sizeof(hdr->nregions)) < 0) goto fail;
    if (hdr->nregions) {
        size_t n = hdr->nregions * sizeof(*hdr->regions);
        hdr->regions = smalloc(n);
        if (hread(fp, hdr->regions, n) < 0) goto fail;
    }

    if (hread(fp, &hdr->region_size, sizeof(hdr->region_size)) < 0) goto fail;
    if (hread(fp, &hdr->nregions_x, sizeof(hdr->nregions_x)) < 0) goto fail;
    if (hread(fp, &hdr->nregions_y, sizeof(hdr->nregions_y)) < 0) goto fail;
    if (Fheader.region_magic[3] == '2') {
        // old filter file
        int n=0;
        if (hread(fp, &n, sizeof(int)) < 0) goto fail;
        hdr->nreads = n;
    } else {
        if (hread(fp, &hdr->nreads, sizeof(hdr->nreads)) < 0) goto fail;
    }
    if (hread(fp, &hdr->readLength, N_READS * sizeof(*hdr->readLength)) < 0) goto fail;

    if (hread(fp, &hdr->filterDataSize, sizeof(hdr->filterDataSize)) < 0) goto fail;
    if (hdr->filterDataSize > 0) {
        hdr->filterData = malloc(hdr->filterDataSize);
        if (hread(fp, hdr->filterData, hdr->filterDataSize) < 0) die("readFilterData() failed\n");
    }

    return hdr;

 fail:
    die("Oops. readHeader() failed\n");
    return NULL;
}

void writeHeader(hFILE *fp, Header *hdr)
{
    if (hwrite(fp, &hdr->lane, sizeof(hdr->lane)) < 0) goto fail;
    if (hwrite(fp, &hdr->coord_shift, sizeof(hdr->coord_shift)) < 0) goto fail;
    if (hwrite(fp, &hdr->coord_factor, sizeof(hdr->coord_factor)) < 0) goto fail;
    if (hwrite(fp, &hdr->ntiles, sizeof(hdr->ntiles)) < 0) goto fail;
    for (int n=0; n < hdr->ntiles; n++) {
        if (hwrite(fp, &hdr->tileArray[n], sizeof(hdr->tileArray[n])) < 0) goto fail;
        if (hwrite(fp, &hdr->tileReadCountArray[n], sizeof(hdr->tileReadCountArray[n])) < 0) goto fail;
    }
    if (hwrite(fp, &hdr->nregions, sizeof(hdr->nregions)) < 0) goto fail;
    if (hwrite(fp, hdr->regions, hdr->nregions * sizeof(*hdr->regions)) < 0) goto fail;
    if (hwrite(fp, &hdr->region_size, sizeof(hdr->region_size)) < 0) goto fail;
    if (hwrite(fp, &hdr->nregions_x, sizeof(hdr->nregions_x)) < 0) goto fail;
    if (hwrite(fp, &hdr->nregions_y, sizeof(hdr->nregions_y)) < 0) goto fail;
    if (hwrite(fp, &hdr->nreads, sizeof(hdr->nreads)) < 0) goto fail;
    if (hwrite(fp, &hdr->readLength, N_READS * sizeof(*hdr->readLength)) < 0) goto fail;
    return;

fail:
    die("writeHeader() failed\n");
}

/*
 * Open and load a list of filter files
 */
static void openFilters(va_t *fnames)
{
    for (int n = 0; n < fnames->end; n++) {
        char *fname = fnames->entries[n];
        hFILE *fp = hopen(fname, "r");
        if (!fp) die("Can't open file %s\n", fname);
        readFheader(fp);
        while (1) {
            Header *hdr = readHeader(fp);
            if (!hdr) break;
            LaneArray[hdr->lane] = hdr;
            hdr->stats_nreads = 0;
            hdr->stats_nfiltered = 0;
        }
        if (hclose(fp)) die("Failed to close filter %s\n", fname);
    }
}

/*
 * Read the supplied SNP (.rod) file into a hash table
 */
static HashTable *readSnpFile(opts_t *opts)
{
    hFILE *fp;
    HashTable *snp_hash;
    static const int line_size = 8192;
    char line[line_size];

    if (!opts->snp_file) return NULL;

    if (opts->verbose) display("reading snp file %s\n", opts->snp_file);

    fp = hopen(opts->snp_file, "rb");
    if (!fp) die("ERROR: can't open known snp file %s: %s\n", opts->snp_file, strerror(errno));

    snp_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS);
    if (snp_hash) die("ERROR: creating snp hash table\n");

    while (hgets(line, line_size, fp)) {
        char key[128];
        HashData hd;
        int bin, start, end;
        char chrom[100];

        if (4 != sscanf(line, "%d\t%99s\t%d\t%d", &bin, chrom, &start, &end)) {
            die("ERROR: reading snp file\n%s\n", line);
        }

        /* N.B rod start is 0 based */
        snprintf(key, sizeof(key), "%s:%d", chrom, start);
        hd.i = 0;
        if (NULL == HashTableAdd(snp_hash, key, strlen(key), hd, NULL)) {
            die("ERROR: building snp hash table\n");
        }
    }

    if (hclose(fp)) die("Can't close SNP file");

    return snp_hash;
}

// create key for region table hash
static void makeRegionKey(region_key_t *key, uint32_t x, uint32_t y)
{
    key->x = x;
    key->y = y;
}

/*
 * initialise a region table entry
 */

static void initialiseRegionTableEntry(RegionTableEntry_t *rt)
{
    rt->align     = 0;
    rt->mismatch  = 0;
    rt->insertion = 0;
    rt->deletion  = 0;
    rt->soft_clip = 0;
    rt->known_snp = 0;
    rt->quality   = 0;
    rt->state     = 0;
}

/*
 * Free the Region Table
 */
static void freeRTS(opts_t *s, Header *hdr, RegionTable_t rts)
{
    int itile, read, cycle;

    for (itile=0; itile<hdr->ntiles; itile++) {
        for (read=0; read<N_READS; read++) {
            if( NULL == rts[itile*N_READS+read]) continue;
            for (cycle=0; cycle < hdr->readLength[read]; cycle++)
                free(rts[itile*N_READS+read][cycle]);
            free(rts[itile*N_READS+read]);
        }
    }
    free(rts);
}

/*
 * Free all the options
 */
static void free_opts(opts_t *opts)
{
    if (!opts) return;
    va_free(opts->filters);
    free(opts->snp_file);
    free(opts->in_bam_file);
    free(opts->working_dir);
    free(opts->output);
    free(opts->apply_stats_out);
    free(opts->tileviz);
    free(opts->argv_list);
    free(opts->input_fmt);
    free(opts->output_fmt);
    free(opts);
}

/*
 * initialise tileviz image
 */
static gdImagePtr initImage(int width, int height, char *base, int type, int read, int cycle, int length)
{
    gdImagePtr im = gdImageCreate(width, height);

    if( NULL == colour_table ){
        colour_table = smalloc(N_COLOURS * sizeof(int));
    }

    gdImageColorAllocate(im, 0, 0, 0); // black - the background colour

    // white + graduated shades of blue from light to dark
    colour_table[COLOUR_LEVEL_0]  = gdImageColorAllocate(im, 255, 255, 255);
    colour_table[COLOUR_LEVEL_1]  = gdImageColorAllocate(im, 211, 222, 235);   
    colour_table[COLOUR_LEVEL_2]  = gdImageColorAllocate(im, 189, 206, 225);
    colour_table[COLOUR_LEVEL_3]  = gdImageColorAllocate(im, 167, 190, 215);
    colour_table[COLOUR_LEVEL_4]  = gdImageColorAllocate(im, 145, 174, 205);
    colour_table[COLOUR_LEVEL_5]  = gdImageColorAllocate(im, 124, 157, 195);
    colour_table[COLOUR_LEVEL_6]  = gdImageColorAllocate(im, 102, 141, 185);
    colour_table[COLOUR_LEVEL_7]  = gdImageColorAllocate(im,  80, 125, 175);
    colour_table[COLOUR_LEVEL_8]  = gdImageColorAllocate(im,  58, 109, 165);
    colour_table[COLOUR_LEVEL_9]  = gdImageColorAllocate(im,  36,  93, 155);
    colour_table[COLOUR_LEVEL_10] = gdImageColorAllocate(im,  15,  77, 146);
    colour_table[COLOUR_LEVEL_11] = gdImageColorAllocate(im,   0,  61, 136);

    // specific colours
    colour_table[COLOUR_TEXT]        = gdImageColorAllocate(im, 239, 239, 239); // light grey
    colour_table[COLOUR_QC_FAIL]     = gdImageColorAllocate(im, 255,   0,   0); // red
    colour_table[COLOUR_ZERO_QUAL]   = gdImageColorAllocate(im, 255,   0,   0); // red
    colour_table[COLOUR_LOW_QUAL]    = gdImageColorAllocate(im, 244, 211,  71); // yellow
    colour_table[COLOUR_MEDIUM_QUAL] = gdImageColorAllocate(im,  21,  58, 144); // dark blue
    colour_table[COLOUR_HIGH_QUAL]   = gdImageColorAllocate(im, 185, 212, 246); // light blue
    
    char str[1024];

    if( NULL != base ) gdImageString(im, gdFontSmall, 3, 1, (unsigned char *)base, colour_table[COLOUR_TEXT]);

    if( cycle < 0 ){
        sprintf(str, "%c_%s", (read == 2 ? 'R' : 'F'), image_names[type]);
    }else{
        sprintf(str, "%0*d%c_%s", length, cycle, (read == 2 ? 'R' : 'F'), image_names[type]);
    }
    gdImageString(im, gdFontSmall, 3, 11, (unsigned char *)str, colour_table[COLOUR_TEXT]);

    return im;
}

/*
 * generate the tileviz report as a HTML file
 */
static void report(opts_t *opts, Header *hdr)
{
    char *base;
    int filename_sz;
    char *filename;
    FILE *fp;
    int image, read, cycle;

    if (0 >= hdr->ntiles)
        return;

    filename_sz = (NULL == opts->tileviz ? 0 : strlen(opts->tileviz)) + 100;
    filename = smalloc(filename_sz);

    sprintf(filename, "%s_lane%d.html", opts->tileviz, hdr->lane);
    fp = fopen(filename, "w+");
    if (!fp) die("Can't open tileviz file %s: %s\n", filename, strerror(errno));

    if (opts->verbose) display("Generating report %s\n", filename);

    base = strrchr(opts->tileviz, '/');
    if( NULL == base )
        base = opts->tileviz;
    else {
        base++;
    }
    
    // initialise the report
    fprintf(fp, "<html>\n");
    fprintf(fp, "<head>\n");
    fprintf(fp, "  <title>Tile Visualisation for %s Lane %d</title>\n", base, hdr->lane);
    fprintf(fp, "  <style type=\"text/css\">\n");
    fprintf(fp, "    table {background-color: rgb(200,200,200)}\n");
    fprintf(fp, "    td {padding: 3px;}\n");
    fprintf(fp, "  </style>\n");
    fprintf(fp, "</head>\n");
    fprintf(fp, "<body>\n");
    fprintf(fp, "  <h3>Tile Visualisation for %s Lane %d</h3>\n", base, hdr->lane);

    // add summary images to the report
    fprintf(fp, "  <h4>Summary</h4>\n");
    fprintf(fp, "  <table>\n");
    fprintf(fp, "    <tr>\n");
    for (image=0; image<N_IMAGES; image++) {
        for (read = 0; read < N_READS; read++) {
            if (0 == hdr->readLength[read]) continue;
            sprintf(filename, "%s_lane%d/%c_%s.png", base, hdr->lane, (read == 2 ? 'R' : 'F'), image_names[image]);
            fprintf(fp, "      <td><img src=\"%s\" /></td>\n", filename);
        }
    }
    fprintf(fp, "    </tr>\n");
    fprintf(fp, "  </table>\n");

    // add cycle by cycle images to the report
    for (read = 0; read < N_READS; read++) {
        if (0 == hdr->readLength[read]) continue;
        int length = (hdr->readLength[read] > 99 ? 3 : (hdr->readLength[read] > 9 ? 2 : 1));
        int image_count = 0;
        fprintf(fp, "  <h4>%s  Read per Cycle</h4>\n", (read == 2 ? "Reverse" : "Forward"));
        fprintf(fp, "  <table>\n");
        fprintf(fp, "    <tr>\n");
        for (cycle = 0; cycle < hdr->readLength[read]; cycle++) {
            for (image=1; image<N_IMAGES; image++) {
	            fprintf(fp, (image ==(N_IMAGES-1) ? "      <td style=\"padding-right:10px;\">" : "      <td>"));
                sprintf(filename, "%s_lane%d/%0*d%c_%s.png", base, hdr->lane, length, cycle+1, (read == 2 ? 'R' : 'F'), image_names[image]);
	            fprintf(fp, "<img src=\"%s\" /></td>\n", filename);
                image_count++;
            }
            // have we reached the end of the row
            if (image_count > NUM_IMAGES_IN_REPORT_ROW) {
                fprintf(fp, "    </tr>\n");
                image_count = 0;
                // are we going to output another row
                if ((cycle+1) < hdr->readLength[read]) {
                    fprintf(fp, "    <tr>\n");
                }
            }
        }
        fprintf(fp, "  </table>\n");
    }

    // finalise the report
    fprintf(fp, "</body>\n");
    fprintf(fp, "</html>\n");

    fclose(fp);
    free(filename);

    return;
}

/*
 * generate tileviz images
*/
static void tileviz(opts_t *opts, Header *hdr, RegionTable_t rts)
{
    int num_surfs = 1;
    int num_cols = 1;
    int num_rows = 1;
    int image_width, image_height;
    gdImagePtr im[N_IMAGES];
    char *base;
    int filename_sz;
    char *filename;
    FILE *fp;
    int image, iregion, ix, iy, read, itile, cycle;

    if (0 >= hdr->ntiles)
        return;

    if (opts->verbose) display("Writing tileviz images to %s_lane%d\n", opts->tileviz, hdr->lane);

    // calculate the number of surfaces, columns and rows, tiles are numbered as follows SCRR where S(surface), C(column) and R(row)
    if( 1 < hdr->ntiles ){
        for (itile=0; itile < hdr->ntiles; itile++) {
            //int tile = s->tileArray[itile];
            int tile = hdr->tileArray[itile];
            int surf = tile / 1000;
            int col = (tile - 1000 * surf) / 100;
            int row = tile % 100;
            num_surfs = max(num_surfs, surf);
            num_cols = max(num_cols, col);
            num_rows = max(num_rows, row);
        }
    }

    image_width = hdr->nregions_x * num_cols * num_surfs + (num_surfs > 1 ? IMAGE_COLUMN_GAP : 0);
    image_height = (hdr->nregions_y + 1) * num_rows + IMAGE_LABEL_HEIGHT;
    
    filename_sz = (NULL == opts->tileviz ? 0 : strlen(opts->tileviz)) + 100;
    filename = smalloc(filename_sz);

    sprintf(filename, "mkdir -p %s_lane%d", opts->tileviz, hdr->lane);
    if (system(filename)) die("Can't make tileviz directory %s_lane%d: %s\n", opts->tileviz, hdr->lane, strerror(errno));

    base = strrchr(opts->tileviz, '/');
    if( NULL == base )
        base = opts->tileviz;
    else {
        base++;
    }
    
    // create the summary images, marking as bad any regions which would be removed or marked as qc failed when the filter is applied
    for (read = 0; read < N_READS; read++) {
        if (0 == hdr->readLength[read]) continue;

        im[IMAGE_COVERAGE]  = initImage(image_width, image_height, base, IMAGE_COVERAGE,  read, -1, 0);
        im[IMAGE_DELETION]  = initImage(image_width, image_height, base, IMAGE_DELETION,  read, -1, 0);
        im[IMAGE_INSERTION] = initImage(image_width, image_height, base, IMAGE_INSERTION, read, -1, 0);
        im[IMAGE_MISMATCH]  = initImage(image_width, image_height, base, IMAGE_MISMATCH,  read, -1, 0);
        im[IMAGE_QUALITY]   = initImage(image_width, image_height, base, IMAGE_QUALITY,   read, -1, 0);

        for (itile=0; itile < hdr->ntiles; itile++) {
            int tile = hdr->tileArray[itile];
            int surf = 1;
            int col = 1;
            int row = 1;

            if( 1 < hdr->ntiles ){
                surf = tile / 1000;
                col = (tile - 1000 * surf) / 100;
                row = tile % 100;
            }

            iregion = 0;
            for (ix = 0; ix < hdr->nregions_x; ix++) {
                for (iy = 0; iy < hdr->nregions_y; iy++) {
                    if (hdr->regions[iregion] >= 0) {
                        RegionTableEntry_t summary_rt;
                        initialiseRegionTableEntry(&summary_rt);

                        int bad_cycle_count = 0;
                        // summary quality is the minimum average quality, initialise to a large value
                        summary_rt.quality = 100.0;
                        for (cycle = 0; cycle < hdr->readLength[read]; cycle++) {
                            RegionTableEntry_t *rt = rts[itile*N_READS+read][cycle] + hdr->regions[iregion];
                            int n = rt->align + rt->insertion + rt->deletion + rt->soft_clip + rt->known_snp;
                            if (0 == n) continue;
                            // coverage should be the same for all cycles
                            summary_rt.align = n;
                            // for quality values calculate an average value
                            rt->quality /= n;
                            // ignore the last cycle of any read which has a higher error rate and lower quality values
                            // ignore first cycle of the reverse read which has a high error rate and lower quality values due to library prep
                            if( (read == 2 && cycle == 0) || (cycle == (hdr->readLength[read]-1)) ) continue;
                            // for mismatch, insertion and deletion take the maximum over all cycles
                            summary_rt.mismatch = max(summary_rt.mismatch, rt->mismatch);
                            summary_rt.insertion = max(summary_rt.insertion, rt->insertion);
                            summary_rt.deletion = max(summary_rt.deletion, rt->deletion);
                            // for quality values take the minimum over all cycles
                            summary_rt.quality = min(summary_rt.quality, rt->quality);
          			        if (rt->state & REGION_STATE_MASK) bad_cycle_count++;
                        }
                        if (bad_cycle_count) summary_rt.state |= REGION_STATE_BAD;

                        int n = summary_rt.align;
                        if (n) {
                            int x = (surf-1) * (hdr->nregions_x * num_cols + IMAGE_COLUMN_GAP) + (col-1) * hdr->nregions_x + ix;
                            int y = IMAGE_LABEL_HEIGHT + (row-1) * (hdr->nregions_y + 1) + iy;
                            int colour = (n > COLOUR_LEVEL_11 ? COLOUR_LEVEL_11 : n);
                            // mark bad regions with COLOUR_QC_FAIL in coverage image
                            if (summary_rt.state & REGION_STATE_BAD) colour = COLOUR_QC_FAIL;
                            gdImageSetPixel(im[IMAGE_COVERAGE],  x, y, colour_table[colour]);
                            // for mismatch, insertion and deletion convert to a percentage and bin 0(<=0), 1(<=10), 2(<=20), ...
                            colour = (10.0 * summary_rt.deletion) / n + (summary_rt.deletion ? 1 : 0);
                            gdImageSetPixel(im[IMAGE_DELETION],  x, y, colour_table[colour]);
                            colour = (10.0 * summary_rt.insertion) / n + (summary_rt.insertion ? 1 : 0);
                            gdImageSetPixel(im[IMAGE_INSERTION], x, y, colour_table[colour]);
                            colour = (10.0 * summary_rt.mismatch) / n + (summary_rt.mismatch ? 1 : 0);
                            gdImageSetPixel(im[IMAGE_MISMATCH],  x, y, colour_table[colour]);
                            // for quality use thresholds >30, >15, >=5 and <5
                            if (summary_rt.quality > 30) {
                                colour = COLOUR_HIGH_QUAL;
                            } else if (summary_rt.quality > 15) {
                                colour = COLOUR_MEDIUM_QUAL;
                            } else if (summary_rt.quality < 5) {
                                colour = COLOUR_ZERO_QUAL;
                            } else {
                                colour = COLOUR_LOW_QUAL;
                            }
                            gdImageSetPixel(im[IMAGE_QUALITY],  x, y, colour_table[colour]);
                        }
                    }

                    iregion++;
                }
            }
        }

        for (image=0; image<N_IMAGES; image++) {
            sprintf(filename, "%s_lane%d/%c_%s.png", opts->tileviz, hdr->lane, (read == 2 ? 'R' : 'F'), image_names[image]);
            fp = fopen(filename, "w+");
            if (!fp) die("Can't open tileviz file %s: %s\n", filename, strerror(errno));
            gdImagePng(im[image], fp);
            fclose(fp);
            gdImageDestroy(im[image]);
        }
    }

    // create cycle by cycle images
    for (read = 0; read < N_READS; read++) {
        if (0 == hdr->readLength[read]) continue;

        int length = (hdr->readLength[read] > 99 ? 3 : (hdr->readLength[read] > 9 ? 2 : 1));
        for (cycle = 0; cycle < hdr->readLength[read]; cycle++) {

            im[IMAGE_DELETION]  = initImage(image_width, image_height, base, IMAGE_DELETION,  read, cycle+1, length);
            im[IMAGE_INSERTION] = initImage(image_width, image_height, base, IMAGE_INSERTION, read, cycle+1, length);
            im[IMAGE_MISMATCH]  = initImage(image_width, image_height, base, IMAGE_MISMATCH,  read, cycle+1, length);
            im[IMAGE_QUALITY]   = initImage(image_width, image_height, base, IMAGE_QUALITY,   read, cycle+1, length);

            for (itile=0; itile < hdr->ntiles; itile++) {
                int tile = hdr->tileArray[itile];
                int surf = 1;
                int col = 1;
                int row = 1;

                if( 1 < hdr->ntiles ){
                    surf = tile / 1000;
                    col = (tile - 1000 * surf) / 100;
                    row = tile % 100;
                }

                iregion = 0;
                for (ix = 0; ix < hdr->nregions_x; ix++) {
                    for (iy = 0; iy < hdr->nregions_y; iy++) {
                        if (hdr->regions[iregion] >= 0) {
                            RegionTableEntry_t *rt = rts[itile*N_READS+read][cycle] + hdr->regions[iregion];
                            int n = rt->align + rt->insertion + rt->deletion + rt->soft_clip + rt->known_snp;
                            if (n) {
                                int x = (surf-1) * (hdr->nregions_x * num_cols + IMAGE_COLUMN_GAP) + (col-1) * hdr->nregions_x + ix;
                                int y = IMAGE_LABEL_HEIGHT + (row-1) * (hdr->nregions_y + 1) + iy;
                                int colour;
                                // for mismatch, insertion and deletion convert to a percentage and bin 0(<=0), 1(<=10), 2(<=20), ...
                                colour = (10.0 * rt->deletion) / n + (rt->deletion ? 1 : 0);
                                gdImageSetPixel(im[IMAGE_DELETION],  x, y, colour_table[colour]);
                                colour = (10.0 * rt->insertion) / n + (rt->insertion ? 1 : 0);
                                gdImageSetPixel(im[IMAGE_INSERTION], x, y, colour_table[colour]);
                                colour = (10.0 * rt->mismatch) / n + (rt->mismatch ? 1 : 0);
                                gdImageSetPixel(im[IMAGE_MISMATCH],  x, y, colour_table[colour]);
                                // for quality use thresholds >30, >15, >=5 and <5
                                if (rt->quality > 30) {
                                    colour = COLOUR_HIGH_QUAL;
                                } else if (rt->quality > 15) {
                                    colour = COLOUR_MEDIUM_QUAL;
                                } else if (rt->quality < 5) {
                                    colour = COLOUR_ZERO_QUAL;
                                } else {
                                    colour = COLOUR_LOW_QUAL;
                                }
                                gdImageSetPixel(im[IMAGE_QUALITY],  x, y, colour_table[colour]);
                            }
                        }
                        iregion++;
                    }
                }
            }
            
            for (image=1; image<N_IMAGES; image++) {
                sprintf(filename, "%s_lane%d/%0*d%c_%s.png", opts->tileviz, hdr->lane, length, cycle+1, (read == 2 ? 'R' : 'F'), image_names[image]);
                fp = fopen(filename, "w+");
                if (!fp) die("Can't open tileviz file %s: %s\n", filename, strerror(errno));
                gdImagePng(im[image], fp);
                fclose(fp);
                gdImageDestroy(im[image]);
            }
        }
    }

    // generate the report
    report(opts, hdr);
    
    return;
}

/*
 * setup a mapping between each potential region and the observed regions
 */
static void regionMapping(Header *hdr)
{
    int iregion, ix, iy;
    region_key_t key;

    if (hdr->nregions <= 0) return;

    hdr->regions = smalloc(hdr->nregions * sizeof(int));
    iregion = 0;
    for (ix = 0; ix < hdr->nregions_x; ix++) {
        for (iy = 0; iy < hdr->nregions_y; iy++) {
            makeRegionKey(&key, ix, iy);
            HashItem *hi = HashTableSearch(hdr->region_hash, (char *)&key, sizeof(region_key_t));
            if (hi) hdr->regions[iregion++] = hi->data.i;
            else    hdr->regions[iregion++] = -1;
        }
    }

    return;
}

/*
 * calc the relative size of the regions we use to set the region state
*/

static int setScaleFactor(opts_t *opts, Header *hdr, RegionTable_t rts)
{
    int scale_factor = 1, region_min_count = 0;
    
	if (0 >= hdr->ntiles)
		return scale_factor;

    // set the region_min_count so that at least 2 reads are required to pass all thresholds
    if ((region_min_count * opts->region_mismatch_threshold) < 2.0 )
        region_min_count = ceil(2.0 / opts->region_mismatch_threshold);
    if ((region_min_count * opts->region_insertion_threshold) < 2.0 )
        region_min_count = ceil(2.0 / opts->region_insertion_threshold);
    if ((region_min_count * opts->region_deletion_threshold) < 2.0 )
        region_min_count = ceil(2.0 / opts->region_deletion_threshold);
    if (opts->verbose) display("State region: region_min_count=%d\n", region_min_count);
    opts->region_min_count = region_min_count;

    int region_size = hdr->region_size;
    int nregions_x = hdr->nregions_x;
    int nregions_y = hdr->nregions_y;
    int nregions = hdr->nregions;

    // what is the average #reads per region, assume coverage is reasonably uniform over the whole lane
    int region_count = (float)hdr->nreads / (float)(hdr->ntiles * nregions);
    if (opts->verbose) display("State region: nregions_x=%d nregions_y=%d nregions=%d region_size=%d region_count=%d\n",
            nregions_x, nregions_y, nregions, region_size, region_count);

    // increase the region size until at the average region count exceeds the minimum region count
    while ( region_count < opts->region_min_count ){
        scale_factor++;
        region_size = scale_factor * hdr->region_size;
        nregions_x = ceil((float)hdr->nregions_x / (float)scale_factor);
        nregions_y = ceil((float)hdr->nregions_y / (float)scale_factor);
        nregions = nregions_x * nregions_y;
        region_count = (float)hdr->nreads / (float)(hdr->ntiles * nregions);
        if (opts->verbose) display("State region: nregions_x=%d nregions_y=%d nregions=%d region_size=%d region_count=%d\n",
                nregions_x, nregions_y, nregions, region_size, region_count);
        // the region size cannot exceed the tile size
        if (nregions == 1) break;
    }

	return scale_factor;
}

/*
 * set the region state
*/

static void setRegionState(opts_t *opts, Header *hdr, RegionTable_t rts)
{
    int scale_factor, nregions_x_state, nregions_y_state, nregions_state;
    RegionTableEntry_t *state_rts = NULL;
	int itile, read, cycle, iregion, ix, iy;

	if (0 >= hdr->ntiles)
		return;

    scale_factor = setScaleFactor(opts, hdr, rts);

    if (scale_factor > 1) {
        if (opts->verbose) display("State region: %dx%d filter regions\n", scale_factor, scale_factor);
        nregions_x_state = ceil((float)hdr->nregions_x / (float)scale_factor);
        nregions_y_state = ceil((float)hdr->nregions_y / (float)scale_factor);
        nregions_state = nregions_x_state * nregions_y_state;
        state_rts = malloc(nregions_state * sizeof(RegionTableEntry_t));
    }else{
        nregions_x_state = hdr->nregions_x;
        nregions_y_state = hdr->nregions_y;
        nregions_state = hdr->nregions;
    }
    
    for (itile=0; itile<hdr->ntiles; itile++) {
        for (read = 0; read < N_READS; read++) {
			for (cycle = 0; cycle < hdr->readLength[read]; cycle++) {
                if (NULL != state_rts) {
                    /* re-initialise the state RT */
                    for (iregion=0; iregion<nregions_state; iregion++)
                        initialiseRegionTableEntry(&state_rts[iregion]);
                    /* fill the state RT */
                    iregion = 0;
                    for (ix = 0; ix < hdr->nregions_x; ix++) {
                        int ix_state = ix / scale_factor;
                        for (iy = 0; iy < hdr->nregions_y; iy++) {
                            int iy_state = iy / scale_factor;
                            int iregion_state = ix_state * nregions_y_state + iy_state;
                            RegionTableEntry_t *state_rt = &state_rts[iregion_state];
                            if (hdr->regions[iregion] >= 0) {
                                RegionTableEntry_t *rt = rts[itile*N_READS+read][cycle] + hdr->regions[iregion];
                                state_rt->align     += rt->align;
                                state_rt->mismatch  += rt->mismatch;
                                state_rt->insertion += rt->insertion;
                                state_rt->deletion  += rt->deletion;
                                state_rt->soft_clip += rt->soft_clip;
                                state_rt->known_snp += rt->known_snp;
                                state_rt->quality   += rt->quality;
                            }
                            iregion++;
                        }
                    }
                }
                /* set the state of the state RT */
                for( iregion = 0; iregion < nregions_state; iregion++) {
                    RegionTableEntry_t *rt;
                    if (NULL != state_rts) {
                        rt = &state_rts[iregion];
                    }else{
                        if (hdr->regions[iregion] < 0) continue;
                        rt = rts[itile*N_READS+read][cycle] + hdr->regions[iregion];
                    }
                    rt->state = 0;
                    // coverage
                    int n = rt->align + rt->insertion + rt->deletion + rt->soft_clip + rt->known_snp;
                    // coverage - mark sparse bins
                    if (n < opts->region_min_count) rt->state |= REGION_STATE_COVERAGE;
                    // correct for sparse bins by assuming ALL bins have atleast region_min_count clusters
                    n = max(n, opts->region_min_count);
                    // mismatch - mark bins with maximum mismatch rate > threshold
                    if (((float)rt->mismatch  / (float)n) >= opts->region_mismatch_threshold)  rt->state |= REGION_STATE_MISMATCH;
                    // insertion - mark bins with maximum insertion rate > threshold
                    if (((float)rt->insertion / (float)n) >= opts->region_insertion_threshold) rt->state |= REGION_STATE_INSERTION;
                    // deletion - mark bins with maximum deletion rate > threshold
                    if (((float)rt->deletion  / (float)n) >= opts->region_deletion_threshold)  rt->state |= REGION_STATE_DELETION;
                }
                if (NULL != state_rts) {
                    /* set the state of the regions using the state RT */
                    iregion = 0;
                    for (ix = 0; ix < hdr->nregions_x; ix++) {
                        int ix_state = ix / scale_factor;
                        for (iy = 0; iy < hdr->nregions_y; iy++) {
                            int iy_state = iy / scale_factor;
                            int iregion_state = ix_state * nregions_y_state + iy_state;
                            RegionTableEntry_t *state_rt = &state_rts[iregion_state];
                            if (hdr->regions[iregion] >= 0) {
                                RegionTableEntry_t *rt = rts[itile*N_READS+read][cycle] + hdr->regions[iregion];
                                rt->state = state_rt->state;
                            }
                            iregion++;
                        }
                    }
                }
			}
		}
	}

    if( NULL != state_rts) free(state_rts);

	// ignoring low coverage, if all regions for each tile/cycle with a non-zero state have the same state
    // and the fraction of regions with this state exceeds a theshold set the state for the whole tile/cycle
    for (itile=0; itile<hdr->ntiles; itile++) {
        for (read = 0; read < N_READS; read++) {
			for (cycle = 0; cycle < hdr->readLength[read]; cycle++) {
				int tile_state = -1, nregions = 0;
                for (iregion=0; iregion<hdr->nregions; iregion++) {
                    if (hdr->regions[iregion] >= 0) {
                        RegionTableEntry_t *rt = rts[itile*N_READS+read][cycle] + hdr->regions[iregion];
                        int state = rt->state & ~REGION_STATE_COVERAGE;
                        if (!state) continue;
                        if (tile_state == -1) tile_state = state;
                        if (state != tile_state) break;
                        nregions++;
                    }
                }
				if (iregion == hdr->nregions && (((float)nregions/(float)hdr->nregions) >= TILE_REGION_THRESHOLD)) {
                    for (iregion=0; iregion<hdr->nregions; iregion++) {
                        if (hdr->regions[iregion] >= 0) {
                            RegionTableEntry_t *rt = rts[itile*N_READS+read][cycle] + hdr->regions[iregion];
                            rt->state = tile_state | (rt->state & REGION_STATE_COVERAGE);
                        }
                    }
                }
			}
		}
	}

    if (!opts->verbose) return;

	// for each tile/cycle output a count of regions with by state
    for (itile=0; itile<hdr->ntiles; itile++) {
		int tile = hdr->tileArray[itile];
        for (read = 0; read < N_READS; read++) {
			for (cycle = 0; cycle < hdr->readLength[read]; cycle++) {
                int mismatch = 0, insertion = 0, deletion = 0, soft_clip = 0;
                long quality_bases = 0, quality_errors = 0;
                for (iregion=0; iregion<hdr->nregions; iregion++) {
                    if (hdr->regions[iregion] >= 0) {
                        RegionTableEntry_t *rt = rts[itile*N_READS+read][cycle] + hdr->regions[iregion];
                        if (rt->state & REGION_STATE_MISMATCH)  mismatch++;
                        if (rt->state & REGION_STATE_INSERTION) insertion++;
                        if (rt->state & REGION_STATE_DELETION)  deletion++;
                        if (rt->state & REGION_STATE_SOFT_CLIP) soft_clip++;
                        quality_bases  += rt->align;
                        quality_errors += rt->mismatch;
                    }
                }
                float ssc = 1.0;
                float quality = -10.0 * log10((quality_errors + ssc)/(quality_bases + ssc));
                if (opts->verbose) 
                    display("tile=%-4d read=%1d cycle=%-3d quality=%.2f mismatch=%-4d insertion=%-4d deletion=%-4d soft_clip=%-4d\n",
                        tile, read, cycle, quality, mismatch, insertion, deletion, soft_clip);
			}
		}
	}
}

/*
* discard the filter if the total number of reads is less than ntiles * MIN_TILE_READ_COUNT
*
* remove individual tiles with less than MIN_TILE_READ_COUNT reads
*/

static void removeBadTiles(Header *hdr)
{
    int ngood_tiles = 0;
    uint64_t nreads = 0, tile_threshold = 0, threshold = 0;
    int itile, read;

	if (0 >= hdr->ntiles) {
        display("No data in filter\n");
		return;
    }
    
    for (itile=0; itile < hdr->ntiles; itile++) {
	    nreads += hdr->tileReadCountArray[itile];
    }

    for (read=0; read < N_READS; read++) {
        if (0 == hdr->readLength[read]) continue;
        tile_threshold += MIN_TILE_READ_COUNT;
    }
    threshold = hdr->ntiles * tile_threshold;
    
	if (0 < nreads && nreads < threshold) {
        display("Discarding filter nreads %" PRIu64 " < %" PRIu64 "\n", nreads, threshold);
		return;
    }

    for (itile=0; itile < hdr->ntiles; itile++) {
        if( hdr->tileReadCountArray[itile] && hdr->tileReadCountArray[itile] < tile_threshold ) {
            display("Discarding filter for tile %d tile_read_count %d < %d\n", hdr->tileArray[itile], hdr->tileReadCountArray[itile], tile_threshold);
            hdr->tileArray[itile] = -1;
        } else {
            ngood_tiles++;
        }
    }

    hdr->ngood_tiles = ngood_tiles;
}

/*
 * Write the filter file to disk
 */
static void writeFilter(opts_t *s, RegionTable_t *rtsArray) 
{
	hFILE *fp;
	int itile, read, cycle, iregion;
	Header *hdr;
    RegionTable_t rts;

    fp = hopen(s->filters->entries[0], "w");
	if (!fp) die("Can't open filter file %s: %s\n", s->filters->entries[0], strerror(errno));

	strncpy(Fheader.region_magic, REGION_MAGIC, sizeof(Fheader.region_magic));
	strncpy(Fheader.cmdLine, s->argv_list, sizeof(Fheader.cmdLine));
    if (hwrite(fp, &Fheader, sizeof(Fheader)) < 0) die("writeFheader() failed\n");;

    for (int lane=1; lane < SF_MAX_LANES; lane++) {
        hdr = LaneArray[lane];
        if (!hdr) continue;
        rts = rtsArray[lane];
        writeHeader(fp,hdr);
            
        int n=0;
        for (int i=0; i < N_READS; i++) n += hdr->readLength[i];
        hdr->filterDataSize = hdr->ntiles * n * hdr->nregions;
        if (hwrite(fp, &hdr->filterDataSize, sizeof(hdr->filterDataSize)) < 0) die("writeHeader() failed\n");;

        for (itile=0; itile<hdr->ntiles; itile++) {
            for (read = 0; read < N_READS; read++) {
                for (cycle = 0; cycle < hdr->readLength[read]; cycle++) {
                    for (iregion=0; iregion<hdr->nregions; iregion++) {
                        int state = 0;
                        if (hdr->regions[iregion] >= 0) {
                            RegionTableEntry_t *rt = rts[itile*N_READS+read][cycle] + hdr->regions[iregion];
                            state = rt->state;
                        }
                        hputc(state, fp);
                    }
                }
            }
        }

    }

	if (hclose(fp)) die("Failed to close \n", s->filters->entries[0]);
}

static int findRegion(opts_t *opts, RegionTable_t rts, Header *hdr, int x, int y)
{
    int iregion = -1;
    uint32_t ix = x2region(hdr, x);
    uint32_t iy = x2region(hdr, y);
    region_key_t key;
    HashItem *hi;

    makeRegionKey(&key, ix, iy);
    hi = HashTableSearch(hdr->region_hash, (char *)&key, sizeof(key));
    if (!hi) {
        HashData hd;
        iregion = hdr->region_hash->nused;
        hd.i = iregion;
        if ( NULL == HashTableAdd(hdr->region_hash, (char *)&key, sizeof(key), hd, NULL) ) {
            die("ERROR: building rts hash table\n");
        }
        int nregions_x = (ix >= hdr->nregions_x ? (ix + 1) : hdr->nregions_x);
        int nregions_y = (iy >= hdr->nregions_y ? (iy + 1) : hdr->nregions_y);
        int nregions = nregions_x * nregions_y;
        if (nregions > hdr->nregions ){
            int itile, read, cycle;
            hdr->nregions_x = nregions_x;
            hdr->nregions_y = nregions_y;
            hdr->nregions = nregions;
            for (itile=0; itile < hdr->ntiles; itile++) {
                for (read=0; read < N_READS; read++) {
                    if (NULL == rts[itile*N_READS+read]) continue;
                    for (cycle=0; cycle < hdr->readLength[read]; cycle++) {
                        int new_iregion;
                        rts[itile*N_READS+read][cycle] = srealloc(rts[itile*N_READS+read][cycle], hdr->nregions * sizeof(RegionTableEntry_t));
                        for (new_iregion = iregion; new_iregion < hdr->nregions; new_iregion++) {
                            RegionTableEntry_t *rt = rts[itile*N_READS+read][cycle] + new_iregion;
                            initialiseRegionTableEntry(rt);
                        }
                    }
                }
            }
        }
        iregion = hd.i;
    }else{
        iregion = hi->data.i;
    }

     return iregion;
}

static void updateRegionTable(Header *hdr, RegionTable_t rts, int read, int iregion, int *read_qual, int *read_mismatch)
{
    /* update region table */
	int cycle;
    for (cycle = 0; cycle < hdr->readLength[read]; cycle++) {
        RegionTableEntry_t *rt = rts[read][cycle] + iregion;
        if (read_mismatch[cycle] & BASE_INSERTION) rt->insertion++;
        if (read_mismatch[cycle] & BASE_DELETION) rt->deletion++;
        if (read_mismatch[cycle] & BASE_SOFT_CLIP) rt->soft_clip++;
        if (read_mismatch[cycle] & BASE_KNOWN_SNP) { 
            rt->known_snp++;
        } else {
            if (read_mismatch[cycle] & BASE_ALIGN) rt->align++;
            if (read_mismatch[cycle] & BASE_MISMATCH) rt->mismatch++;
        }
        rt->quality += read_qual[cycle];
    }

	return;
}

/*
 * create an ordered array of tiles and re-order the RegionTable by tile
 */
static RegionTable_t orderRegionTableByTile(opts_t *s, Header *hdr, RegionTable_t rts)
{
    RegionTable_t new_rts = NULL;
	int itile, read;

    if (0 >= hdr->ntiles)
        return new_rts;

    // create a sorted array of tiles
	int *tileArray = smalloc(hdr->ntiles * sizeof(int));
	size_t *tileReadCountArray = smalloc(hdr->ntiles * sizeof(size_t));
	for (itile=0; itile < hdr->ntiles; itile++)
        tileArray[itile] = hdr->tileArray[itile];
    qsort(tileArray, hdr->ntiles, sizeof(int), int_sort);

	// re-order the region table by tile
    new_rts = smalloc(hdr->ntiles * N_READS * sizeof(RegionTableEntry_t **));
	for (itile=0; itile < hdr->ntiles; itile++) {
	    int tile = tileArray[itile];
        int old_itile = tile2index(hdr, tile);
   	    for (read=0; read < N_READS; read++) {
            tileReadCountArray[itile] = hdr->tileReadCountArray[old_itile];
            new_rts[itile*N_READS+read] = rts[old_itile*N_READS+read];
        }
    }

    free(rts);
    free(hdr->tileArray); hdr->tileArray = tileArray;
    free(hdr->tileReadCountArray); hdr->tileReadCountArray = tileReadCountArray;
    return new_rts;
}

/*
 * Takes the bam file as input and updates the region table
 *
 * Assumption: within a single input file, all reads are the same length and
 * we're using unclipped data.
 *
 * Returns: 0 written for success
 *	   -1 for failure
 */
static RegionTable_t *makeRegionTable(opts_t *opts, BAMit_t *fp_bam)
{
    RegionTable_t rts = NULL;
    RegionTable_t *rtsArray = smalloc(SF_MAX_LANES * sizeof(rts));
    Header *hdr;

	static const int bam_read_buff_size = 1024;
	char bam_read_seq[bam_read_buff_size];
	int bam_read_qual[bam_read_buff_size];
	int bam_read_mismatch[bam_read_buff_size];

	bam1_t *bam;

    for (int lane=0; lane < SF_MAX_LANES; lane++) rtsArray[lane] = NULL;

	/* loop over reads in the bam file */
	while (1) {
		int bam_lane = -1, bam_tile = -1, bam_read = -1, bam_x = -1, bam_y = -1, read_length;

		bam = parse_bam_readinfo(fp_bam, &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_read, NULL);
        if (!bam) {
			break;	/* break on end of BAM file */
		}

		if (BAM_FUNMAP & bam->core.flag) continue;
		if (BAM_FQCFAIL & bam->core.flag) continue;
		if (BAM_FSECONDARY & bam->core.flag) continue;
		if (BAM_FSUPPLEMENTARY & bam->core.flag) continue;
		if (BAM_FPAIRED & bam->core.flag) {
			if (0 == (BAM_FPROPER_PAIR & bam->core.flag)) {
				continue;
			}
		}
        
        rts = rtsArray[bam_lane];
        hdr = LaneArray[bam_lane];
        if (!hdr) {
            hdr = sf_initHdr();
            hdr->region_size = opts->region_size;
            hdr->lane = bam_lane;
            LaneArray[bam_lane] = hdr;
        }
        read_length = bam->core.l_qseq;
        if (0 == hdr->readLength[bam_read]) {
            hdr->readLength[bam_read] = read_length;
        }

        if (hdr->readLength[bam_read] != read_length) {
            fprintf(stderr,
                    "Error: inconsistent read lengths "
                    "within bam file for read %d.\n"
                    "have length %ld, previously it was %d.\n",
                    bam_read, (long) read_length, hdr->readLength[bam_read]);
            exit(EXIT_FAILURE);
        }

		parse_bam_alignments(fp_bam, bam, bam_read_seq, bam_read_qual, NULL, bam_read_mismatch,
                                                  bam_read_buff_size, opts->snp_hash);

        // lookup tile in tile array
        int itile = tile2index(hdr,bam_tile);
	    if (itile < 0) {
            itile = hdr->ntiles;
            hdr->ntiles++;
            hdr->tileArray = srealloc(hdr->tileArray, hdr->ntiles * sizeof(int));
            hdr->tileReadCountArray = srealloc(hdr->tileReadCountArray, hdr->ntiles * sizeof(size_t));
            hdr->tileArray[itile] = bam_tile;
            hdr->tileReadCountArray[itile] = 0;
            rts = srealloc(rts, hdr->ntiles * N_READS * sizeof(RegionTableEntry_t **));
            rtsArray[bam_lane] = rts;
            for (int read=0;read<N_READS;read++)
                rts[itile*N_READS+read] = NULL;
            if (opts->verbose) fprintf(stderr, "Processing lane %d tile %i (%" PRIu64 ")\n", bam_lane, bam_tile, hdr->nreads);
        }
        hdr->tileReadCountArray[itile]++;

        if (NULL == rts[itile*N_READS+bam_read]) {
            int cycle, iregion;
            rts[itile*N_READS+bam_read] = smalloc(read_length * sizeof(RegionTableEntry_t *));
            for(cycle=0;cycle<read_length;cycle++) {
                rts[itile*N_READS+bam_read][cycle] = smalloc(hdr->nregions * sizeof(RegionTableEntry_t));
                for(iregion=0;iregion<hdr->nregions;iregion++) {
                    RegionTableEntry_t *rt = rts[itile*N_READS+bam_read][cycle] + iregion;
                    initialiseRegionTableEntry(rt);
                }
            }
        }

        int iregion = findRegion(opts, rts, hdr, bam_x, bam_y);
        updateRegionTable(hdr, &rts[itile*N_READS], bam_read, iregion, bam_read_qual, bam_read_mismatch);

        hdr->nreads++;
	}

	bam_destroy1(bam);

    /* re-order the RegionTable by tile */
    for (int lane=1; lane < SF_MAX_LANES; lane++) {
        Header *hdr = LaneArray[lane];
        if (!hdr) continue;
	    rtsArray[lane] = orderRegionTableByTile(opts, hdr, rtsArray[lane]);
    }

    /* setup a mapping between each potential region and the observed regions */
    for (int lane=1; lane < SF_MAX_LANES; lane++) {
        if (LaneArray[lane]) regionMapping(LaneArray[lane]);
    }

    return rtsArray;
}

/*
 * Takes a bam file as input and outputs a filtered bam file
 *
 * Returns: 0 written for success
 *	   -1 for failure
 */
static int filter_bam(opts_t *opts, BAMit_t *fp_in_bam, BAMit_t *fp_out_bam)
{
	bam1_t *bam;
    bool ignore = false;

	/* loop over reads in the input bam file */
	while (1) {
        /* apply filter */

        int bam_lane = -1, bam_tile = -1, bam_read = -1, bam_x = -1, bam_y = -1;
        ignore = false;

        bam = parse_bam_readinfo(fp_in_bam, &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_read, NULL);
        if (!bam) {
            break;	/* break on end of BAM file */
        }

        Header *hdr = LaneArray[bam_lane];
        if (hdr) {
            if (hdr->ngood_tiles) {
                int iregion = xy2region(hdr, bam_x, bam_y);
                char *state = getFilterData(hdr, tile2index(hdr,bam_tile), 0, 0, iregion);
                if (state != NULL) {
                    int read, cycle, bad_cycle_count = 0;
                    for (read = 0; read < N_READS; read++) {
                        int maxCycle = hdr->readLength[read];
                        for (cycle = 0; cycle < maxCycle; cycle++) {
                            if (*state & REGION_STATE_MASK)
                                bad_cycle_count++;
                            state += hdr->nregions;
                        }
                    }

                    if (bad_cycle_count) {
                        hdr->stats_nfiltered++;
                        if (opts->qcfail) 
                            bam->core.flag |= BAM_FQCFAIL;
                        else
                            ignore = true;
                    }
                }
            }
            hdr->stats_nreads++;
        }

		if (!ignore) {
            if (0 > sam_write1(fp_out_bam->f, fp_out_bam->h, bam)) die("Error: writing bam file\n");
        }
	}

	bam_destroy1(bam);

	return 0;
}


static void usage(FILE *usagefp)
{
	fprintf(usagefp, "Usage: bambi spatial_filter [command] [options] bam_file\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "Command: must be one and only one of:\n");
	fprintf(usagefp, " -D               Dump filter file in ascii text format to stdout\n");
	fprintf(usagefp, " -c               Create filter file from BAM file\n");
	fprintf(usagefp, " -a               Apply filter file to a BAM file\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "Options:\n");
	fprintf(usagefp, " -v --verbose     display progress messages to stderr\n");
	fprintf(usagefp, "    --output-fmt  BAM output format [sam|bam|cram] [default: bam]\n");
	fprintf(usagefp, "    --input-fmt   BAM input format [sam|bam|cram] [default: bam]\n");
	fprintf(usagefp, "    --compression-level\n");
    fprintf(usagefp, "                  Compression level for output BAM\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "Comand specific options:\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "    all other commands require:\n");
	fprintf(usagefp, "      -F --filter file\n");
	fprintf(usagefp, "                  Filter filename e.g. 8088.filter\n");
	fprintf(usagefp, "                  no default: must be supplied\n");
	fprintf(usagefp, "                  or\n");
	fprintf(usagefp, "                  comma separated list of filter files (for apply filter only).\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "    create filter:\n");
	fprintf(usagefp, "      -s --snp_file file\n");
	fprintf(usagefp, "                 set of snps to be removed\n");
	fprintf(usagefp, "                 file in Reference Ordered Data (ROD) format\n");
	fprintf(usagefp, "      --region_size\n");
	fprintf(usagefp, "                 default %d\n", REGION_SIZE);
	fprintf(usagefp, "      --region_mismatch_threshold\n");
	fprintf(usagefp, "                 threshold for setting region mismatch state\n");
	fprintf(usagefp, "                 default %-6.4f\n", REGION_MISMATCH_THRESHOLD);
	fprintf(usagefp, "      --region_insertion_threshold\n");
	fprintf(usagefp, "                 threshold for setting region insertion state\n");
	fprintf(usagefp, "                 default %-6.4f\n", REGION_INSERTION_THRESHOLD);
	fprintf(usagefp, "      --region_deletion_threshold\n");
	fprintf(usagefp, "                 threshold for setting region deletion state\n");
	fprintf(usagefp, "                 default %-6.4f\n", REGION_DELETION_THRESHOLD);
	fprintf(usagefp, "      -t prefix\n");
	fprintf(usagefp, "                 generate tileviz files in this directory\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "    apply filter:\n");
	fprintf(usagefp, "      -o         output\n");
	fprintf(usagefp, "                 Output bam file name\n");
	fprintf(usagefp, "                 default: stdout\n");
	fprintf(usagefp, "      -f         mark filtered reads as QCFAIL\n");
	fprintf(usagefp, "                 default: do not output filtered reads\n");
	fprintf(usagefp, "      -l         apply_stats\n");
	fprintf(usagefp, "                 apply status message output\n");
	fprintf(usagefp, "                 default: stderr\n");
	fprintf(usagefp, "\n");
}

/*
 * Create filter command
 */
static void calculateFilter(opts_t *opts)
{
	BAMit_t *fp_input_bam;
	RegionTable_t *rtsArray;
    
	fp_input_bam = BAMit_open(opts->in_bam_file, 'r', opts->input_fmt, 0, NULL);
	if (NULL == fp_input_bam) {
		die("ERROR: can't open bam file %s: %s\n", opts->in_bam_file, strerror(errno));
	}

    /* read the snp_file */
    opts->snp_hash = readSnpFile(opts);

	rtsArray = makeRegionTable(opts, fp_input_bam);

	/* close the bam file */
	BAMit_free(fp_input_bam);

	if (opts->verbose) {
        uint64_t traces = 0;
        for (int n=1; n < SF_MAX_LANES; n++) if (LaneArray[n]) traces += LaneArray[n]->nreads;
		display("Processed %" PRIu64 " traces\n", traces);
		if (opts->snp_hash) {
			size_t nsnps = 0;
			int ibucket;
			for (ibucket = 0; ibucket < opts->snp_hash->nbuckets; ibucket++) {
				HashItem *hi;
				for (hi = opts->snp_hash->bucket[ibucket]; hi; hi = hi->next)
					if (hi->data.i) nsnps += hi->data.i;
			}
			display("Ignored %lu snps\n", nsnps);
		}
	}

    for (int lane=1; lane < SF_MAX_LANES; lane++) {
        if (LaneArray[lane]) setRegionState(opts, LaneArray[lane], rtsArray[lane]);
    }

    if (!opts->filters) {
        display("Writing filter to stdout\n");
        va_push(opts->filters,"/dev/stdout");
    }
    writeFilter(opts, rtsArray);

    if (opts->tileviz) {
        for (int n=1; n < SF_MAX_LANES; n++) {
            if (LaneArray[n]) tileviz(opts, LaneArray[n], rtsArray[n]);
        }
    }
    
    for (int lane=1; lane < SF_MAX_LANES; lane++) {
	    if (LaneArray[lane]) freeRTS(opts, LaneArray[lane], rtsArray[lane]);
    }
    free(rtsArray);
}

/*
 * Apply filter command
 */
static void applyFilter(opts_t *s)
{
	BAMit_t *fp_input_bam;
	BAMit_t *fp_output_bam;
	hFILE *apply_stats_fd = NULL;
	char *out_bam_file = NULL;
	char *apply_stats_file = NULL;

    openFilters(s->filters);

    /* remove bad tiles from region table */
    for (int lane=1; lane < SF_MAX_LANES; lane++) {
        Header *hdr = LaneArray[lane];
        if (hdr) removeBadTiles(hdr);
    }

    // Create output BAM filename
    out_bam_file = smalloc(strlen(s->working_dir) + strlen(s->output) + 16);
    out_bam_file[0] = 0;
    if (s->output[0] != '/') {
        strcpy(out_bam_file, s->working_dir);
        strcat(out_bam_file, "/");
    }
    strcat(out_bam_file, s->output);

    // Create stats filename
    apply_stats_file = smalloc(strlen(s->working_dir) + strlen(s->apply_stats_out) + 16);
    apply_stats_file[0] = 0;
    if (s->apply_stats_out[0] != '/') {
        strcpy(apply_stats_file, s->working_dir);
        strcat(apply_stats_file, "/");
    }
    strcat(apply_stats_file, s->apply_stats_out);

	fp_input_bam = BAMit_open(s->in_bam_file, 'r', s->input_fmt, 0, NULL);
	if (NULL == fp_input_bam) {
		die("ERROR: can't open bam file %s: %s\n", s->in_bam_file, strerror(errno));
	}


	fp_output_bam = BAMit_open(out_bam_file, 'w', s->output_fmt, s->compression_level, NULL);
	if (NULL == fp_output_bam) {
		die("ERROR: can't open bam file %s: %s\n", out_bam_file, strerror(errno));
	}
    // copy input to output header
    bam_hdr_destroy(fp_output_bam->h); fp_output_bam->h = bam_hdr_dup(fp_input_bam->h);

	char concat_cmd[2048];
	strcat(concat_cmd, s->argv_list);
	bam_header_add_pg("spf", "spatial_filter", "A program to apply a spatial filter", concat_cmd, fp_output_bam->h);
    if (sam_hdr_write(fp_output_bam->f, fp_output_bam->h) < 0) die("Can't write %s header\n", out_bam_file);
	free(out_bam_file);


	if (-1 == filter_bam(s, fp_input_bam, fp_output_bam)) {
		die("ERROR: failed to filter bam file %s\n", s->in_bam_file);
	}

	BAMit_free(fp_input_bam);
	BAMit_free(fp_output_bam);

//
//  Display stats
//
	if (NULL == (apply_stats_fd=hopen(apply_stats_file, "w"))) {
		die("ERROR: failed to open apply status log %s\n", apply_stats_file);
	}

    for (int n=0; n < SF_MAX_LANES; n++) {
        Header *hdr = LaneArray[n];
        if (!hdr) continue;
        char buffer[64];
        sprintf(buffer, "Lane %d\t", hdr->lane);
        hputs(buffer, apply_stats_fd);
        sprintf(buffer, "Processed %" PRIu64 " \t", hdr->stats_nreads);
        hputs(buffer, apply_stats_fd);
        sprintf(buffer, "%s %" PRIu64 " traces\n", (s->qcfail ? "Failed" : "Removed"), hdr->stats_nfiltered);
        hputs(buffer, apply_stats_fd);
    }
    hputs("\n", apply_stats_fd);

	if (hclose(apply_stats_fd)) die("Can't close stats file");
    free(apply_stats_file);
}

/*
 * Dump filter file command
 */
static void dumpFilterFile(opts_t *opts)
{
    openFilters(opts->filters);

	printf("Magic:          %s\n", Fheader.region_magic);
	printf("Command Line:   %s\n", Fheader.cmdLine);

    for (int laneNo=1; laneNo < SF_MAX_LANES; laneNo++) {
        Header *hdr = LaneArray[laneNo];;
        if (!hdr) continue;

        printf("\n");
        printf("Lane:           %-5d\n", hdr->lane);
        printf("Coord Shift:    %-5d\n", hdr->coord_shift);
        printf("Coord Factor:   %-5d\n", hdr->coord_shift);
        printf("Region Size:    %-5d\n", hdr->region_size);
        printf("Num Regions:    %-5d\n", hdr->nregions);
        printf("Num Regions X:  %-5d\n", hdr->nregions_x);
        printf("Num Regions Y:  %-5d\n", hdr->nregions_y);
        for (int n=0; n < hdr->nregions; n++) printf("%d ", hdr->regions[n]); 
        printf("\n");
        printf("Num Tiles:      %-5d\n", (int)hdr->ntiles);
        for (int i=0; i < hdr->ntiles; i++) printf("%-5d %-12lu", hdr->tileArray[i], hdr->tileReadCountArray[i]);
        printf("\n");
        printf("Read Length:    ");
        for (int i=0; i < N_READS; i++) printf("%-5d ", hdr->readLength[i]);
        printf("\n");
        printf("Filter Size:    %" PRIu32 "\n", hdr->filterDataSize);

        if (opts->verbose) {
            int itile, read, cycle, region;
            for (itile=0; itile<hdr->ntiles; itile++) {
                for (read=0; read < N_READS; read++) {
                    for (cycle=0; cycle<hdr->readLength[read]; cycle++) {
                        char* state = getFilterData(hdr, itile, read, cycle, 0);
                        if (state != NULL) {
                            for (region=0; region<hdr->nregions; region++) {
                                if (*state & REGION_STATE_MASK)
                                    printf("filtering tile=%d read=%d cycle=%d region=%d\n", hdr->tileArray[itile], read, cycle, region);
                                state++;
                            }
                        }
                    }
                }
            }
        }
    }

}

va_t *parseList(char *arg)
{
    va_t *va = va_init(5,free);
    parse_tags(va,arg);
    return va;
}

/*
 * Parse the command line arguments into a form we can use
 */
opts_t* spatial_filter_parse_args(int argc, char *argv[])
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char *optstring = "vdcafuDF:b:e:o:l:i:p:s:r:R:x:y:t:z:qh?";

	static const struct option lopts[] = {
        {"snp_file", 1, 0, 's'},
        {"snp-file", 1, 0, 's'},
        {"help", 0, 0, 'h'},
        {"filter", 1, 0, 'F'},
        {"rg", 1, 0, 'R'},
        {"verbose", 0, 0, 'v'},
        {"region-size", 1, 0, 'r'},
        {"region_mismatch_threshold", 1, 0, 'z'},
        {"region_insertion_threshold", 1, 0, 'b'},
        {"region_deletion_threshold", 1, 0, 'e'},
        {"output-fmt", 1, 0, 0},
        {"input-fmt", 1, 0, 0},
        {"compression-level", 1, 0, 0},
        {"tileviz", 1, 0, 't'},
        {0, 0, 0, 0}
    };

    opts_t *opts = calloc(sizeof(opts_t), 1);
    if (!opts) { perror("cannot allocate option parsing memory"); return NULL; }

    opts->argv_list = stringify_argv(argc+1, argv-1);
    if (opts->argv_list[strlen(opts->argv_list)-1] == ' ') opts->argv_list[strlen(opts->argv_list)-1] = 0;

    // Set default values
	opts->region_size = REGION_SIZE;
	opts->region_mismatch_threshold = REGION_MISMATCH_THRESHOLD;
	opts->region_insertion_threshold = REGION_INSERTION_THRESHOLD;
	opts->region_deletion_threshold = REGION_DELETION_THRESHOLD;
    opts->filters = NULL;

    int opt;
    int ncmd = 0;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, optstring, lopts, &option_index)) != -1) {
        const char *arg;
        switch (opt) {
            case 'R': break;
			case 'D': opts->dumpFilter = 1; ncmd++;       break;
			case 'c': opts->calculate = 1; ncmd++;        break;
            case 't': opts->tileviz = strdup(optarg);     break;
			case 'a': opts->apply = 1; ncmd++;            break;
			case 'f': opts->qcfail = 1;		              break;
			case 'o': opts->output = strdup(optarg);      break;
			case 's': opts->snp_file = strdup(optarg);    break;
			case 'F': opts->filters = parseList(optarg);  break;
			case 'r': opts->region_size = atoi(optarg);   break;
			case 'z': opts->region_mismatch_threshold = atof(optarg); break;
			case 'b': opts->region_insertion_threshold = atof(optarg); break;
			case 'e': opts->region_deletion_threshold = atof(optarg); break;
			case 'v': opts->verbose = 1;			        break;
			case 'l': opts->apply_stats_out = strdup(optarg);	break;
			case 'h': usage(stdout); free_opts(opts); return NULL;
            case 0:   arg = lopts[option_index].name;
                          if (strcmp(arg, "output-fmt") == 0)              opts->output_fmt = strdup(optarg);
                     else if (strcmp(arg, "input-fmt") == 0)               opts->input_fmt = strdup(optarg);
                     else if (strcmp(arg, "compression-level") == 0)       opts->compression_level = *optarg;
                     else {
                         fprintf(stderr,"\nUnknown option: %s\n\n", arg);
                         usage(stderr); free_opts(opts);
                         return NULL;
                     }
                     break;
            default: fprintf(stderr,"Unknown option: '%c'\n", opt);
                     usage(stderr); free_opts(opts);
                     return NULL;

        }
    }

    if (ncmd > 1) {
        fprintf(stderr,"ERROR: More than one command given\n");
        usage(stderr); free_opts(opts); return NULL; 
    }

    if (ncmd == 0) {
        fprintf(stderr,"ERROR: No command given\n");
        usage(stderr); free_opts(opts); return NULL; 
    }

	if (optind < argc) opts->in_bam_file = strdup(argv[optind]);

	if (!opts->in_bam_file && !opts->dumpFilter) die("Error: no BAM file specified\n");

    if (!opts->filters && (opts->dumpFilter || opts->apply)) die("Error: no filter file specified\n");

	if (opts->calculate) {
   	    if (opts->region_size < 1) die("Error: invalid region size\n");
    }

    if (!opts->calculate && opts->tileviz) display("Warning: no tileviz images will be produced\n");

    if (!opts->apply_stats_out) opts->apply_stats_out = strdup("/dev/stderr");
    if (!opts->output) opts->output = strdup("/dev/stdout");

    return opts;
}

static int spatial_filter(opts_t *opts)
{
	/* preserve starting directory */
	opts->working_dir = getcwd(NULL,0);
	if (!opts->working_dir) {
		die("ERROR: can't obtain working directory: %s\n", strerror(errno));
	}

	/* Dump the filter file */
	if (opts->dumpFilter) dumpFilterFile(opts);

	/* calculate the filter */
    if (opts->calculate) calculateFilter(opts);

	/* apply the  filter */
	if (opts->apply) applyFilter(opts);

	return EXIT_SUCCESS;

}

/*
 * Called from bambi to perform spatial filtering
 *
 * parse the command line arguments, then call the main process
 *
 * returns 0 on success, 1 if there was a problem
 */
int main_spatial_filter(int argc, char *argv[])
{
    int ret = 1;
    opts_t *opts = spatial_filter_parse_args(argc, argv);
    if (opts) ret = spatial_filter(opts);
    free_opts(opts);
    return ret;
}

