/*  spatial_filter.c - This code looks for spatial features given an aligned bam file

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
#include "bamit.h"
#include "hash_table.h"
#include "rts.h"
#include "parse_bam.h"
#include "array.h"
#include "parse.h"

#define PHRED_QUAL_OFFSET  33  // phred quality values offset

#define REGION_MISMATCH_THRESHOLD   0.016  // threshold for setting region mismatch state
#define REGION_INSERTION_THRESHOLD  0.016  // threshold for setting region insertion state
#define REGION_DELETION_THRESHOLD   0.016  // threshold for setting region deletion state

#define TILE_REGION_THRESHOLD  0.75  // threshold for setting region state at tile level

#define MIN_TILE_READ_COUNT  1000 // min number of aligned reads on a tile

#define REGION_STATE_MASK  (REGION_STATE_INSERTION | REGION_STATE_DELETION)  // region mask used to filter reads

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

int *colour_table = NULL;

typedef struct {
	va_t *filters;
    va_t *rgids;
	char *snp_file;
	char *in_bam_file;
	HashTable *snp_hash;
	int *tileArray;
	size_t *tileReadCountArray;
	char *working_dir;
	char *output;
	char *apply_stats_out;
	int read_length[3];
	int calculate;
    bool dumpFilter;
	char *tileviz;
	int apply;
	int qcfail;
	int verbose;
	HashTable *region_hash;
	int region_size;
	int nregions_x;
	int nregions_y;
	int nregions;
	int *regions;
	int region_min_count;
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
 * Read the supplied SNP (.rod) file into a hash table
 */
static HashTable *readSnpFile(opts_t *opts)
{
    FILE *fp;
    HashTable *snp_hash;
    static const int line_size = 8192;
    char line[line_size];

    if (!opts->snp_file) return NULL;

    if (opts->verbose) display("reading snp file %s\n", opts->snp_file);

    fp = fopen(opts->snp_file, "rb");
    if (!fp) die("ERROR: can't open known snp file %s: %s\n", opts->snp_file, strerror(errno));

    snp_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS);
    if (snp_hash) die("ERROR: creating snp hash table\n");

    while (fgets(line, line_size, fp)) {
        char key[100];
        HashData hd;
        int bin, start, end;
        char chrom[100];

        if (4 != sscanf(line, "%d\t%s\t%d\t%d", &bin, chrom, &start, &end)) {
            die("ERROR: reading snp file\n%s\n", line);
        }

        /* N.B rod start is 0 based */
        snprintf(key, sizeof(key), "%s:%d", chrom, start);
        hd.i = 0;
        if (NULL == HashTableAdd(snp_hash, key, strlen(key), hd, NULL)) {
            die("ERROR: building snp hash table\n");
        }
    }

    fclose(fp);

    return snp_hash;
}

/*
 * Convert an integer to a string
 */
static char *_itoa(char *cp, int i)
{
    int j;

    if (i == 0) {
        *cp++ = '0';
        return cp;
    }

    if (i < 0) {
        *cp++ = '-';
        if (i == INT_MIN) {
            *cp++ = '2'; *cp++ = '1'; *cp++ = '4'; *cp++ = '7';
            *cp++ = '4'; *cp++ = '8'; *cp++ = '3'; *cp++ = '6';
            *cp++ = '4'; *cp++ = '8';
            return cp;
        }
        i = -i;
    }

    //if (i < 10)         goto b0;
    if (i < 100)        goto b1;
    //if (i < 1000)       goto b2;
    if (i < 10000)      goto b3;
    //if (i < 100000)     goto b4;
    if (i < 1000000)    goto b5;
    //if (i < 10000000)   goto b6;
    if (i < 100000000)  goto b7;

     if ((j = i / 1000000000)) {*cp++ = j + '0'; i -= j*1000000000; goto x8;}
     if ((j = i / 100000000))  {*cp++ = j + '0'; i -= j*100000000;  goto x7;}
 b7: if ((j = i / 10000000))   {*cp++ = j + '0'; i -= j*10000000;   goto x6;}
     if ((j = i / 1000000))    {*cp++ = j + '0', i -= j*1000000;    goto x5;}
 b5: if ((j = i / 100000))     {*cp++ = j + '0', i -= j*100000;     goto x4;}
     if ((j = i / 10000))      {*cp++ = j + '0', i -= j*10000;      goto x3;}
 b3: if ((j = i / 1000))       {*cp++ = j + '0', i -= j*1000;       goto x2;}
     if ((j = i / 100))        {*cp++ = j + '0', i -= j*100;        goto x1;}
 b1: if ((j = i / 10))         {*cp++ = j + '0', i -= j*10;         goto x0;}
     if (i)                     *cp++ = i + '0';
    return cp;

 x8: *cp++ = i / 100000000 + '0', i %= 100000000;
 x7: *cp++ = i / 10000000  + '0', i %= 10000000;
 x6: *cp++ = i / 1000000   + '0', i %= 1000000;
 x5: *cp++ = i / 100000    + '0', i %= 100000;
 x4: *cp++ = i / 10000     + '0', i %= 10000;
 x3: *cp++ = i / 1000      + '0', i %= 1000;
 x2: *cp++ = i / 100       + '0', i %= 100;
 x1: *cp++ = i / 10        + '0', i %= 10;
 x0: *cp++ = i             + '0';

    return cp;
}

static void makeRegionKey(char *key, int x, int y)
{
    key = _itoa(key, x);
    *key++ = ':';
    key = _itoa(key, y);
    *key = 0;
}

/*
 * initialise a region table entry
 */

static void initialiseRegionTable(RegionTable *rt)
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
static void freeRTS(opts_t *s, int ntiles, RegionTable ***rts)
{
    int itile, read, cycle;

    for (itile=0; itile<ntiles; itile++) {
        for (read=0; read<N_READS; read++) {
            if( NULL == rts[itile*N_READS+read]) continue;
            for (cycle=0; cycle < s->read_length[read]; cycle++)
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
    va_free(opts->rgids);
    free(opts->snp_file);
    free(opts->in_bam_file);
    free(opts->tileArray);
    free(opts->tileReadCountArray);
    free(opts->regions);
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
static void report(opts_t *opts, int ntiles)
{
    char *base;
    int filename_sz;
    char *filename;
    FILE *fp;
    int image, read, cycle;

    if (0 >= ntiles)
        return;

    filename_sz = (NULL == opts->tileviz ? 0 : strlen(opts->tileviz)) + 100;
    filename = smalloc(filename_sz);

    sprintf(filename, "%s.html", opts->tileviz);
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
    fprintf(fp, "  <title>Tile Visualisation for %s</title>\n", base);
    fprintf(fp, "  <style type=\"text/css\">\n");
    fprintf(fp, "    table {background-color: rgb(200,200,200)}\n");
    fprintf(fp, "    td {padding: 3px;}\n");
    fprintf(fp, "  </style>\n");
    fprintf(fp, "</head>\n");
    fprintf(fp, "<body>\n");
    fprintf(fp, "  <h3>Tile Visualisation for %s</h3>\n", base);

    // add summary images to the report
    fprintf(fp, "  <h4>Summary</h4>\n");
    fprintf(fp, "  <table>\n");
    fprintf(fp, "    <tr>\n");
    for (image=0; image<N_IMAGES; image++) {
        for (read = 0; read < N_READS; read++) {
            if (0 == opts->read_length[read]) continue;
            sprintf(filename, "%s/%c_%s.png", base, (read == 2 ? 'R' : 'F'), image_names[image]);
            fprintf(fp, "      <td><img src=\"%s\" /></td>\n", filename);
        }
    }
    fprintf(fp, "    </tr>\n");
    fprintf(fp, "  </table>\n");

    // add cycle by cycle images to the report
    for (read = 0; read < N_READS; read++) {
        if (0 == opts->read_length[read]) continue;
        int length = (opts->read_length[read] > 99 ? 3 : (opts->read_length[read] > 9 ? 2 : 1));
        int image_count = 0;
        fprintf(fp, "  <h4>%s  Read per Cycle</h4>\n", (read == 2 ? "Reverse" : "Forward"));
        fprintf(fp, "  <table>\n");
        fprintf(fp, "    <tr>\n");
        for (cycle = 0; cycle < opts->read_length[read]; cycle++) {
            for (image=1; image<N_IMAGES; image++) {
	            fprintf(fp, (image ==(N_IMAGES-1) ? "      <td style=\"padding-right:10px;\">" : "      <td>"));
                sprintf(filename, "%s/%0*d%c_%s.png", base, length, cycle+1, (read == 2 ? 'R' : 'F'), image_names[image]);
	            fprintf(fp, "<img src=\"%s\" /></td>\n", filename);
                image_count++;
            }
            // have we reached the end of the row
            if (image_count > NUM_IMAGES_IN_REPORT_ROW) {
                fprintf(fp, "    </tr>\n");
                image_count = 0;
                // are we going to output another row
                if ((cycle+1) < opts->read_length[read]) {
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

static void tileviz(opts_t *s, int ntiles, RegionTable ***rts)
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

    if (0 >= ntiles)
        return;

    if (s->verbose) display("Writing tileviz images to %s\n", s->tileviz);

    // calculate the number of surfaces, columns and rows, tiles are numbered as follows SCRR where S(surface), C(column) and R(row)
    if( 1 < ntiles ){
        for (itile=0; itile<ntiles; itile++) {
            int tile = s->tileArray[itile];
            int surf = tile / 1000;
            int col = (tile - 1000 * surf) / 100;
            int row = tile % 100;
            num_surfs = max(num_surfs, surf);
            num_cols = max(num_cols, col);
            num_rows = max(num_rows, row);
        }
    }

    image_width = s->nregions_x * num_cols * num_surfs + (num_surfs > 1 ? IMAGE_COLUMN_GAP : 0);
    image_height = (s->nregions_y + 1) * num_rows + IMAGE_LABEL_HEIGHT;
    
    filename_sz = (NULL == s->tileviz ? 0 : strlen(s->tileviz)) + 100;
    filename = smalloc(filename_sz);

    sprintf(filename, "mkdir -p %s", s->tileviz);
    if (system(filename)) die("Can't make tileviz directory %s: %s\n", s->tileviz, strerror(errno));

    base = strrchr(s->tileviz, '/');
    if( NULL == base )
        base = s->tileviz;
    else {
        base++;
    }
    
    // create the summary images, marking as bad any regions which would be removed or marked as qc failed when the filter is applied
    for (read = 0; read < N_READS; read++) {
        if (0 == s->read_length[read]) continue;

        im[IMAGE_COVERAGE]  = initImage(image_width, image_height, base, IMAGE_COVERAGE,  read, -1, 0);
        im[IMAGE_DELETION]  = initImage(image_width, image_height, base, IMAGE_DELETION,  read, -1, 0);
        im[IMAGE_INSERTION] = initImage(image_width, image_height, base, IMAGE_INSERTION, read, -1, 0);
        im[IMAGE_MISMATCH]  = initImage(image_width, image_height, base, IMAGE_MISMATCH,  read, -1, 0);
        im[IMAGE_QUALITY]   = initImage(image_width, image_height, base, IMAGE_QUALITY,   read, -1, 0);

        for (itile=0; itile<ntiles; itile++) {
            int tile = s->tileArray[itile];
            int surf = 1;
            int col = 1;
            int row = 1;

            if( 1 < ntiles ){
                surf = tile / 1000;
                col = (tile - 1000 * surf) / 100;
                row = tile % 100;
            }

            iregion = 0;
            for (ix = 0; ix < s->nregions_x; ix++) {
                for (iy = 0; iy < s->nregions_y; iy++) {
                    if (s->regions[iregion] >= 0) {
                        RegionTable summary_rt;
                        initialiseRegionTable(&summary_rt);

                        int bad_cycle_count = 0;
                        // summary quality is the minimum average quality, initialise to a large value
                        summary_rt.quality = 100.0;
                        for (cycle = 0; cycle < s->read_length[read]; cycle++) {
                            RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                            int n = rt->align + rt->insertion + rt->deletion + rt->soft_clip + rt->known_snp;
                            if (0 == n) continue;
                            // coverage should be the same for all cycles
                            summary_rt.align = n;
                            // for quality values calculate an average value
                            rt->quality /= n;
                            // ignore the last cycle of any read which has a higher error rate and lower quality values
                            // ignore first cycle of the reverse read which has a high error rate and lower quality values due to library prep
                            if( (read == 2 && cycle == 0) || (cycle == (s->read_length[read]-1)) ) continue;
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
                            int x = (surf-1) * (s->nregions_x * num_cols + IMAGE_COLUMN_GAP) + (col-1) * s->nregions_x + ix;
                            int y = IMAGE_LABEL_HEIGHT + (row-1) * (s->nregions_y + 1) + iy;
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
            sprintf(filename, "%s/%c_%s.png", s->tileviz, (read == 2 ? 'R' : 'F'), image_names[image]);
            fp = fopen(filename, "w+");
            if (!fp) die("Can't open tileviz file %s: %s\n", filename, strerror(errno));
            gdImagePng(im[image], fp);
            fclose(fp);
            gdImageDestroy(im[image]);
        }
    }

    // create cycle by cycle images
    for (read = 0; read < N_READS; read++) {
        if (0 == s->read_length[read]) continue;

        int length = (s->read_length[read] > 99 ? 3 : (s->read_length[read] > 9 ? 2 : 1));
        for (cycle = 0; cycle < s->read_length[read]; cycle++) {

            im[IMAGE_DELETION]  = initImage(image_width, image_height, base, IMAGE_DELETION,  read, cycle+1, length);
            im[IMAGE_INSERTION] = initImage(image_width, image_height, base, IMAGE_INSERTION, read, cycle+1, length);
            im[IMAGE_MISMATCH]  = initImage(image_width, image_height, base, IMAGE_MISMATCH,  read, cycle+1, length);
            im[IMAGE_QUALITY]   = initImage(image_width, image_height, base, IMAGE_QUALITY,   read, cycle+1, length);

            for (itile=0; itile<ntiles; itile++) {
                int tile = s->tileArray[itile];
                int surf = 1;
                int col = 1;
                int row = 1;

                if( 1 < ntiles ){
                    surf = tile / 1000;
                    col = (tile - 1000 * surf) / 100;
                    row = tile % 100;
                }

                iregion = 0;
                for (ix = 0; ix < s->nregions_x; ix++) {
                    for (iy = 0; iy < s->nregions_y; iy++) {
                        if (s->regions[iregion] >= 0) {
                            RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                            int n = rt->align + rt->insertion + rt->deletion + rt->soft_clip + rt->known_snp;
                            if (n) {
                                int x = (surf-1) * (s->nregions_x * num_cols + IMAGE_COLUMN_GAP) + (col-1) * s->nregions_x + ix;
                                int y = IMAGE_LABEL_HEIGHT + (row-1) * (s->nregions_y + 1) + iy;
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
                sprintf(filename, "%s/%0*d%c_%s.png", s->tileviz, length, cycle+1, (read == 2 ? 'R' : 'F'), image_names[image]);
                fp = fopen(filename, "w+");
                if (!fp) die("Can't open tileviz file %s: %s\n", filename, strerror(errno));
                gdImagePng(im[image], fp);
                fclose(fp);
                gdImageDestroy(im[image]);
            }
        }
    }

    // generate the report
    report(s, ntiles);
    
    return;
}

/*
 * setup a mapping between each potential region and the observed regions
 */
static void regionMapping(opts_t *s)
{
    int iregion, ix, iy;
    char key[100];

    free(s->regions);
    s->regions = NULL;
    if (s->nregions <= 0) return;

    s->regions = smalloc(s->nregions * sizeof(int));
    iregion = 0;
    for (ix = 0; ix < s->nregions_x; ix++) {
        for (iy = 0; iy < s->nregions_y; iy++) {
            makeRegionKey(key, ix, iy);
            HashItem *hi = HashTableSearch(s->region_hash, key, strlen(key));
            if (hi) s->regions[iregion++] = hi->data.i;
            else    s->regions[iregion++] = -1;
        }
    }

    return;
}

/*
 * calc the relative size of the regions we use to set the region state
*/

static int setScaleFactor(opts_t *s, int ntiles, size_t nreads, RegionTable ***rts)
{
    int scale_factor = 1, region_min_count = 0;
    
	if (0 >= ntiles)
		return scale_factor;

    // set the region_min_count so that at least 2 reads are required to pass all thresholds
    if ((region_min_count * s->region_mismatch_threshold) < 2.0 )
        region_min_count = ceil(2.0 / s->region_mismatch_threshold);
    if ((region_min_count * s->region_insertion_threshold) < 2.0 )
        region_min_count = ceil(2.0 / s->region_insertion_threshold);
    if ((region_min_count * s->region_deletion_threshold) < 2.0 )
        region_min_count = ceil(2.0 / s->region_deletion_threshold);
    if (s->verbose) display("State region: region_min_count=%d\n", region_min_count);
    s->region_min_count = region_min_count;

    int region_size = s->region_size;
    int nregions_x = s->nregions_x;
    int nregions_y = s->nregions_y;
    int nregions = s->nregions;

    // what is the average #reads per region, assume coverage is reasonably uniform over the whole lane
    int region_count = (float)nreads / (float)(ntiles * nregions);
    if (s->verbose) display("State region: nregions_x=%d nregions_y=%d nregions=%d region_size=%d region_count=%d\n",
            nregions_x, nregions_y, nregions, region_size, region_count);

    // increase the region size until at the average region count exceeds the minimum region count
    while ( region_count < s->region_min_count ){
        scale_factor++;
        region_size = scale_factor * s->region_size;
        nregions_x = ceil((float)s->nregions_x / (float)scale_factor);
        nregions_y = ceil((float)s->nregions_y / (float)scale_factor);
        nregions = nregions_x * nregions_y;
        region_count = (float)nreads / (float)(ntiles * nregions);
        if (s->verbose) display("State region: nregions_x=%d nregions_y=%d nregions=%d region_size=%d region_count=%d\n",
                nregions_x, nregions_y, nregions, region_size, region_count);
        // the region size cannot exceed the tile size
        if (nregions == 1) break;
    }

	return scale_factor;
}

/*
 * set the region state
*/

static void setRegionState(opts_t *s, int ntiles, size_t nreads, RegionTable ***rts)
{
    int scale_factor, nregions_x_state, nregions_y_state, nregions_state;
    RegionTable *state_rts = NULL;
	int itile, read, cycle, iregion, ix, iy;

	if (0 >= ntiles)
		return;

    scale_factor = setScaleFactor(s, ntiles, nreads, rts);

    if (scale_factor > 1) {
        if (s->verbose) display("State region: %dx%d filter regions\n", scale_factor, scale_factor);
        nregions_x_state = ceil((float)s->nregions_x / (float)scale_factor);
        nregions_y_state = ceil((float)s->nregions_y / (float)scale_factor);
        nregions_state = nregions_x_state * nregions_y_state;
        state_rts = malloc(nregions_state * sizeof(RegionTable));
    }else{
        nregions_x_state = s->nregions_x;
        nregions_y_state = s->nregions_y;
        nregions_state = s->nregions;
    }
    
    for (itile=0; itile<ntiles; itile++) {
        for (read = 0; read < N_READS; read++) {
			for (cycle = 0; cycle < s->read_length[read]; cycle++) {
                if (NULL != state_rts) {
                    /* re-initialise the state RT */
                    for (iregion=0; iregion<nregions_state; iregion++)
                        initialiseRegionTable(&state_rts[iregion]);
                    /* fill the state RT */
                    iregion = 0;
                    for (ix = 0; ix < s->nregions_x; ix++) {
                        int ix_state = ix / scale_factor;
                        for (iy = 0; iy < s->nregions_y; iy++) {
                            int iy_state = iy / scale_factor;
                            int iregion_state = ix_state * nregions_y_state + iy_state;
                            RegionTable *state_rt = &state_rts[iregion_state];
                            if (s->regions[iregion] >= 0) {
                                RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
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
                    RegionTable *rt;
                    if (NULL != state_rts) {
                        rt = &state_rts[iregion];
                    }else{
                        if (s->regions[iregion] < 0) continue;
                        rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                    }
                    rt->state = 0;
                    // coverage
                    int n = rt->align + rt->insertion + rt->deletion + rt->soft_clip + rt->known_snp;
                    // coverage - mark sparse bins
                    if (n < s->region_min_count) rt->state |= REGION_STATE_COVERAGE;
                    // correct for sparse bins by assuming ALL bins have atleast region_min_count clusters
                    n = max(n, s->region_min_count);
                    // mismatch - mark bins with maximum mismatch rate > threshold
                    if (((float)rt->mismatch  / (float)n) >= s->region_mismatch_threshold)  rt->state |= REGION_STATE_MISMATCH;
                    // insertion - mark bins with maximum insertion rate > threshold
                    if (((float)rt->insertion / (float)n) >= s->region_insertion_threshold) rt->state |= REGION_STATE_INSERTION;
                    // deletion - mark bins with maximum deletion rate > threshold
                    if (((float)rt->deletion  / (float)n) >= s->region_deletion_threshold)  rt->state |= REGION_STATE_DELETION;
                }
                if (NULL != state_rts) {
                    /* set the state of the regions using the state RT */
                    iregion = 0;
                    for (ix = 0; ix < s->nregions_x; ix++) {
                        int ix_state = ix / scale_factor;
                        for (iy = 0; iy < s->nregions_y; iy++) {
                            int iy_state = iy / scale_factor;
                            int iregion_state = ix_state * nregions_y_state + iy_state;
                            RegionTable *state_rt = &state_rts[iregion_state];
                            if (s->regions[iregion] >= 0) {
                                RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
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
    for (itile=0; itile<ntiles; itile++) {
        for (read = 0; read < N_READS; read++) {
			for (cycle = 0; cycle < s->read_length[read]; cycle++) {
				int tile_state = -1, nregions = 0;
                for (iregion=0; iregion<s->nregions; iregion++) {
                    if (s->regions[iregion] >= 0) {
                        RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                        int state = rt->state & ~REGION_STATE_COVERAGE;
                        if (!state) continue;
                        if (tile_state == -1) tile_state = state;
                        if (state != tile_state) break;
                        nregions++;
                    }
                }
				if (iregion == s->nregions && (((float)nregions/(float)s->nregions) >= TILE_REGION_THRESHOLD)) {
                    for (iregion=0; iregion<s->nregions; iregion++) {
                        if (s->regions[iregion] >= 0) {
                            RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                            rt->state = tile_state | (rt->state & REGION_STATE_COVERAGE);
                        }
                    }
                }
			}
		}
	}

    if (!s->verbose) return;

	// for each tile/cycle output a count of regions with by state
    for (itile=0; itile<ntiles; itile++) {
		int tile = s->tileArray[itile];
        for (read = 0; read < N_READS; read++) {
			for (cycle = 0; cycle < s->read_length[read]; cycle++) {
                int mismatch = 0, insertion = 0, deletion = 0, soft_clip = 0;
                long quality_bases = 0, quality_errors = 0;
                for (iregion=0; iregion<s->nregions; iregion++) {
                    if (s->regions[iregion] >= 0) {
                        RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
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
                if (s->verbose) 
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
    size_t nreads = 0, tile_threshold = 0, threshold = 0;
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
        display("Discarding filter nreads %lu < %lu\n", nreads, threshold);
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
static void printFilter(opts_t *s, int ntiles, RegionTable ***rts) 
{
	FILE *fp;
	int itile, read, cycle, iregion;
	Header hdr;

    fp = fopen(s->filters->entries[0], "w+");
	if (!fp) die("Can't open filter file %s: %s\n", s->filters->entries[0], strerror(errno));

	hdr.region_magic = strdup(REGION_MAGIC);
	hdr.coord_shift = COORD_SHIFT; 
	hdr.coord_factor = COORD_FACTOR;
	hdr.ntiles = ntiles;
	if( ntiles > 0 ) {
        hdr.tileArray = s->tileArray;
        hdr.tileReadCountArray = s->tileReadCountArray;
    }
    hdr.nreads = N_READS;
	hdr.region_size = s->region_size;
	hdr.nregions = s->nregions;
	hdr.nregions_x = s->nregions_x;
	hdr.nregions_y = s->nregions_y;
	for (read=0; read < hdr.nreads; read++)
		hdr.readLength[read] = s->read_length[read];
	hdr.cmdLine = strdup(s->argv_list);
	hdr.ncomments = 0;
	writeHeader(fp,&hdr);
        
	free(hdr.region_magic);
	free(hdr.cmdLine);

    for (itile=0; itile<ntiles; itile++) {
		for (read = 0; read < N_READS; read++) {
			for (cycle = 0; cycle < s->read_length[read]; cycle++) {
                for (iregion=0; iregion<s->nregions; iregion++) {
                    int state = 0;
                    if (s->regions[iregion] >= 0) {
                        RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                        state = rt->state;
                    }
					fputc(state, fp);
                }
			}
		}
	}

	fclose(fp);
}

static int findRegion(opts_t *s, RegionTable ***rts, int ntiles, int x, int y)
{
    int iregion = -1;
    int ix = x2region(x, s->region_size);
    int iy = x2region(y, s->region_size);
    char key[100];
    char *cp;
    HashItem *hi;

    makeRegionKey(key, ix, iy);
    hi = HashTableSearch(s->region_hash, key, strlen(key));
    if (!hi) {
        HashData hd;
        iregion = s->region_hash->nused;
        hd.i = iregion;
        if ( NULL == HashTableAdd(s->region_hash, key, strlen(key), hd, NULL) ) {
            die("ERROR: building rts hash table\n");
        }
        int nregions_x = (ix >= s->nregions_x ? (ix + 1) : s->nregions_x);
        int nregions_y = (iy >= s->nregions_y ? (iy + 1) : s->nregions_y);
        int nregions = nregions_x * nregions_y;
        if (nregions > s->nregions ){
            int itile, read, cycle;
            s->nregions_x = nregions_x;
            s->nregions_y = nregions_y;
            s->nregions = nregions;
            for (itile=0; itile < ntiles; itile++) {
                for (read=0; read < N_READS; read++) {
                    if (NULL == rts[itile*N_READS+read]) continue;
                    for (cycle=0; cycle < s->read_length[read]; cycle++) {
                        int new_iregion;
                        rts[itile*N_READS+read][cycle] = srealloc(rts[itile*N_READS+read][cycle], s->nregions * sizeof(RegionTable));
                        for (new_iregion = iregion; new_iregion < s->nregions; new_iregion++) {
                            RegionTable *rt = rts[itile*N_READS+read][cycle] + new_iregion;
                            initialiseRegionTable(rt);
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

static void updateRegionTable(opts_t *s, RegionTable ***rts, int read, int iregion, int *read_qual, int *read_mismatch)
{
    /* update region table */
	int cycle;
    for (cycle = 0; cycle < s->read_length[read]; cycle++) {
        RegionTable *rt = rts[read][cycle] + iregion;
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
static RegionTable ***orderRegionTableByTile(opts_t *s, int *tiles, size_t *nreads, int ntiles, RegionTable ***rts)
{
    RegionTable ***new_rts = NULL;
	int itile, read;

    if (0 >= ntiles)
        return new_rts;

    // create a sorted array of tiles
	s->tileArray = smalloc(ntiles * sizeof(int));
	s->tileReadCountArray = smalloc(ntiles * sizeof(size_t));
	for (itile=0; itile < ntiles; itile++)
        s->tileArray[itile] = tiles[itile];
    qsort(s->tileArray, ntiles, sizeof(int), int_sort);

	// re-order the region table by tile
    new_rts = smalloc(ntiles * N_READS * sizeof(RegionTable **));
	for (itile=0; itile < ntiles; itile++) {
	    int tile = s->tileArray[itile];
        size_t nelem = ntiles;
        void *pitile = lfind(&tile, tiles, &nelem, sizeof(int), &int_cmp);
        int old_itile = ((int*)pitile - tiles);
   	    for (read=0; read < N_READS; read++) {
            s->tileReadCountArray[itile] = nreads[old_itile];
            new_rts[itile*N_READS+read] = rts[old_itile*N_READS+read];
        }
    }

    free(rts);

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
static RegionTable ***makeRegionTable(opts_t *s, BAMit_t *fp_bam, int *bam_ntiles, size_t *bam_nreads)
{
    RegionTable ***rts = NULL;

    int *tiles = NULL;
    size_t *tileReadCounts = NULL;
    int ntiles = 0;

	size_t nreads = 0;

	int lane = -1;

	static const int bam_read_buff_size = 1024;
	char bam_read_seq[bam_read_buff_size];
	int bam_read_qual[bam_read_buff_size];
	int bam_read_mismatch[bam_read_buff_size];

	bam1_t *bam;

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
        
        read_length = bam->core.l_qseq;
        if (0 == s->read_length[bam_read]) {
            s->read_length[bam_read] = read_length;
        }

        if (s->read_length[bam_read] != read_length) {
            fprintf(stderr,
                    "Error: inconsistent read lengths "
                    "within bam file for read %d.\n"
                    "have length %ld, previously it was %d.\n",
                    bam_read, (long) read_length, s->read_length[bam_read]);
            exit(EXIT_FAILURE);
        }

        if (lane == -1) {
            lane = bam_lane;
        }
        if (bam_lane != lane){
            die("Error: Inconsistent lane: have %d previously it was %d.\n", bam_lane, lane);
        }

		parse_bam_alignments(fp_bam, bam, bam_read_seq, bam_read_qual, NULL, bam_read_mismatch,
                                                  bam_read_buff_size, s->snp_hash);

        // lookup tile in tile array
        size_t nelem = ntiles;
	    void *pitile;
	    int itile;
        pitile = lfind(&bam_tile, tiles, &nelem, sizeof(int), &int_cmp);
	    if (NULL == pitile) {
            int read;
            itile = ntiles;
            ntiles++;
            tiles = srealloc(tiles, ntiles * sizeof(int));
            tileReadCounts = srealloc(tileReadCounts, ntiles * sizeof(size_t));
            tiles[itile] = bam_tile;
            tileReadCounts[itile] = 0;
            rts = srealloc(rts, ntiles * N_READS * sizeof(RegionTable **));
            for (read=0;read<N_READS;read++)
                rts[itile*N_READS+read] = NULL;
            if (s->verbose) fprintf(stderr, "Processing tile %i (%lu)\n", bam_tile, nreads);
        }else{
            itile = ((int*)pitile - tiles);
        }
        tileReadCounts[itile]++;

        if (NULL == rts[itile*N_READS+bam_read]) {
            int cycle, iregion;
            rts[itile*N_READS+bam_read] = smalloc(read_length * sizeof(RegionTable *));
            for(cycle=0;cycle<read_length;cycle++) {
                rts[itile*N_READS+bam_read][cycle] = smalloc(s->nregions * sizeof(RegionTable));
                for(iregion=0;iregion<s->nregions;iregion++) {
                    RegionTable *rt = rts[itile*N_READS+bam_read][cycle] + iregion;
                    initialiseRegionTable(rt);
                }
            }
        }

        int iregion = findRegion(s, rts, ntiles, bam_x, bam_y);
        updateRegionTable(s, &rts[itile*N_READS], bam_read, iregion, bam_read_qual, bam_read_mismatch);

        nreads++;
	}

	bam_destroy1(bam);

    /* re-order the RegionTable by tile */
	rts = orderRegionTableByTile(s, tiles, tileReadCounts, ntiles, rts);

    /* setup a mapping between each potential region and the observed regions */
    regionMapping(s);

    if (NULL != tiles) free(tiles);
    if (NULL != tileReadCounts) free(tileReadCounts);

    *bam_ntiles = ntiles;
	*bam_nreads = nreads;

    return rts;
}

/*
 * Takes a bam file as input and outputs a filtered bam file
 *
 * Returns: 0 written for success
 *	   -1 for failure
 */
static int filter_bam(opts_t *s, BAMit_t *fp_in_bam, BAMit_t *fp_out_bam)
{
	int lane = -1;
    char *rgid;
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

        // select filter file depending on RG tag
        if ( (s->rgids->end == 1) && (strcmp(s->rgids->entries[0],"null")==0)) {
            rgid = "null";
        } else {
            uint8_t *p = bam_aux_get(bam,"RG");
            rgid = p ? bam_aux2Z(p) : "null";
        }

        if (setCurrentHdr(rgid)) {
            if (getHdrngood_tiles()) {
                int iregion = xy2region(bam_x, bam_y);
                char* state = getFilterData(bam_tile, 0, 0, iregion);
                if (state != NULL) {
                    int read, cycle, bad_cycle_count = 0;
                    for (read = 0; read < N_READS; read++) {
                        int maxCycle = getHdrReadLength(read);
                        for (cycle = 0; cycle < maxCycle; cycle++) {
                            if (*state & REGION_STATE_MASK)
                                bad_cycle_count++;
                            state += getHdrnregions();
                        }
                    }

                    if (bad_cycle_count) {
                        incHdrStatsnfiltered();
                        if (s->qcfail) 
                            bam->core.flag |= BAM_FQCFAIL;
                        else
                            ignore = true;
                    }
                }
            }
            incHdrStatsnreads();
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
    fprintf(usagefp, "                  There must be a corresponding list of RG IDs (see --rg)\n\n");
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
    fprintf(usagefp, "      -R --rg    Comma separated list of RG IDs. Must correspond to list of filter files\n");
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


static void calculateFilter(opts_t *opts)
{
	BAMit_t *fp_input_bam;
	int ntiles = 0;
	size_t nreads = 0;

	RegionTable ***rts = NULL;
    
	fp_input_bam = BAMit_open(opts->in_bam_file, 'r', opts->input_fmt, 0);
	if (NULL == fp_input_bam) {
		die("ERROR: can't open bam file %s: %s\n", opts->in_bam_file, strerror(errno));
	}

    /* read the snp_file */
    opts->snp_hash = readSnpFile(opts);

    opts->region_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE|HASH_FUNC_JENKINS);

	rts = makeRegionTable(opts, fp_input_bam, &ntiles, &nreads);

	/* close the bam file */
	BAMit_free(fp_input_bam);

	if (opts->verbose) {
		display("Processed %8lu traces\n", nreads);
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

	/* back to where we belong */
	if (chdir(opts->working_dir)) {
        fprintf(stderr,"chdir(%s) failed\n", opts->working_dir);
        exit(EXIT_FAILURE);
    }

    setRegionState(opts, ntiles, nreads, rts);

    if (!opts->filters) {
        display("Writing filter to stdout\n");
        va_push(opts->filters,"/dev/stdout");
    }
    printFilter(opts, ntiles, rts);

    if (opts->tileviz) tileviz(opts, ntiles, rts);
    
    HashTableDestroy(opts->region_hash, 0);
	freeRTS(opts, ntiles, rts);
}

static void applyFilter(opts_t *s)
{
	BAMit_t *fp_input_bam;
	BAMit_t *fp_output_bam;
	FILE *apply_stats_fd = NULL;
	char out_mode[5] = "wb";
	char *out_bam_file = NULL;
	char *apply_stats_file = NULL;
	FILE *fp;
    int read;

    openFilters(s->filters,s->rgids);

    /* remove bad tiles from region table */
    for (int n=0; n < s->rgids->end; n++) {
        char *rgid = s->rgids->entries[n];
        Header *hdr = getHdr(rgid);
        if (!hdr) die("Can't find filter file for %s\n", rgid);
        removeBadTiles(hdr);
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

	fp_input_bam = BAMit_open(s->in_bam_file, 'r', s->input_fmt, 0);
	if (NULL == fp_input_bam) {
		die("ERROR: can't open bam file %s: %s\n", s->in_bam_file, strerror(errno));
	}


	fp_output_bam = BAMit_open(out_bam_file, 'w', s->output_fmt, s->compression_level);
	if (NULL == fp_output_bam) {
		die("ERROR: can't open bam file %s: %s\n", out_bam_file, strerror(errno));
	}
    // copy input to output header
    bam_hdr_destroy(fp_output_bam->h); fp_output_bam->h = bam_hdr_dup(fp_input_bam->h);

	char concat_cmd[2048];
    for (int n=0; n < s->rgids->end; n++) {
        char *rgid = s->rgids->entries[n];
        Header *hdr = getHdr(rgid);
	    strcpy(concat_cmd, hdr->cmdLine);
	    strcat(concat_cmd, " ; ");
    }
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
	if (NULL == (apply_stats_fd=fopen(apply_stats_file, "w"))) {
		die("ERROR: failed to open apply status log %s\n", apply_stats_file);
	}

    va_t *va = HdrHash2Array();
    for (int n=0; n < va->end; n++) {
        Header *hdr = (Header *)va->entries[n];
        fprintf(apply_stats_fd, "%s\t", (strcmp(hdr->rgid,"null") ? hdr->rgid : "Total"));
        fprintf(apply_stats_fd, "Processed %" PRIu64 " \t", hdr->stats_nreads);
        fprintf(apply_stats_fd, "%s %" PRIu64 " traces\n", (s->qcfail ? "Failed" : "Removed"), hdr->stats_nfiltered);
    }
    fprintf(apply_stats_fd, "\n");

	fclose(apply_stats_fd);

	/* back to where we belong */
	if (chdir(s->working_dir)) {
        fprintf(stderr,"chdir(%s) failed\n", s->working_dir);
        exit(EXIT_FAILURE);
    }
}

static void dumpFilterFile(opts_t *opts)
{
	Header *hdr;
	int i;

    openFilters(opts->filters, opts->rgids);
	hdr = getHdr(NULL);
    setCurrentHdr(NULL);

	printf("Magic:          %s\n", hdr->region_magic);
	printf("Coord Shift:    %-5d\n", hdr->coord_shift);
	printf("Coord Factor:   %-5d\n", hdr->coord_shift);
	printf("Region Size:    %-5d\n", hdr->region_size);
	printf("Num Regions:    %-5d\n", hdr->nregions);
	printf("Num Tiles:      %-5d\n", hdr->ntiles);
	for (i=0; i < hdr->ntiles; i++) printf("%-5d %-12lu", hdr->tileArray[i], hdr->tileReadCountArray[i]);
	printf("\n");
	printf("Read Length:    ");
	for (i=0; i < hdr->nreads; i++) printf("%-5d ", hdr->readLength[i]);
	printf("\n");
	printf("Command Line:   %s\n", hdr->cmdLine);
	for (i=0; i < hdr->ncomments; i++) {
		if (i) printf("                %s\n", hdr->comments[i]);
		else   printf("Comments:       %s\n", hdr->comments[i]);
	}

	if (opts->verbose) {
        int tile, read, cycle, region;
        for (tile=0; tile<hdr->ntiles; tile++)
            for (read=0; read<hdr->nreads; read++)
                for (cycle=0; cycle<hdr->readLength[read]; cycle++) {
                    char* state = getFilterData(hdr->tileArray[tile], read, cycle, 0);
                    if (state != NULL) {
                        for (region=0; region<hdr->nregions; region++) {
                            if (*state & REGION_STATE_MASK)
                                printf("filtering tile=%d read=%d cycle=%d region=%d\n", hdr->tileArray[tile], read, cycle, region);
                            state++;
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
    opts->rgids = NULL;

    int opt;
    int ncmd = 0;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, optstring, lopts, &option_index)) != -1) {
        const char *arg;
        switch (opt) {
			case 'D': opts->dumpFilter = 1; ncmd++;       break;
			case 'c': opts->calculate = 1; ncmd++;        break;
            case 't': opts->tileviz = strdup(optarg);     break;
			case 'a': opts->apply = 1; ncmd++;            break;
			case 'f': opts->qcfail = 1;		              break;
			case 'o': opts->output = strdup(optarg);      break;
			case 's': opts->snp_file = strdup(optarg);    break;
			case 'F': opts->filters = parseList(optarg);  break;
            case 'R': opts->rgids = parseList(optarg);    break;
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

    if (!opts->rgids) opts->rgids = parseList("null");
    if (opts->rgids->end != opts->filters->end) die("Not the same number of filter files and RG IDs\n");

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

