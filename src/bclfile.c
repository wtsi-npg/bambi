/* bclfile.c

   Functions to read and parse an Illumina BCL file.

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
#include <libgen.h>

#include "bclfile.h"

#define BCL_BASE_ARRAY "ACGT"
#define BCL_UNKNOWN_BASE 'N'

// return most significant digit of tile
int bcl_tile2surface(int tile)
{
    int surface = tile/1000;
    while (surface > 9) surface = surface/10;
    return surface;
}

static int uncompressBlock(char* abSrc, int nLenSrc, char* abDst, int nLenDst )
{
    z_stream zInfo ={0};
    zInfo.total_in=  zInfo.avail_in=  nLenSrc;
    zInfo.total_out= zInfo.avail_out= nLenDst;
    zInfo.next_in= abSrc;
    zInfo.next_out= abDst;

    int nErr, nRet= -1;
    nErr= inflateInit2( &zInfo, 15+32 );               // zlib function
    if (nErr != Z_OK) fprintf(stderr,"inflateInit() failed: %d\n", nErr);
    if ( nErr == Z_OK ) {
        nErr= inflate( &zInfo, Z_FINISH );     // zlib function
        if ( (nErr == Z_STREAM_END) || (nErr == Z_BUF_ERROR) ) {
            nRet= zInfo.total_out;
            if (nRet != nLenDst) fprintf(stderr,"inflate() returned %d: expected %d\n",nErr,nLenDst);
            nErr = Z_OK;
        }
        if (nErr != Z_OK) {
            fprintf(stderr,"inflate() returned: %d\n", nErr);
            fprintf(stderr,"avail_in=%d  avail_out=%d  total_out=%ld\n", zInfo.avail_in, zInfo.avail_out, zInfo.total_out);
        }
    }
    inflateEnd( &zInfo );   // zlib function
    return( nRet ); // -1 or len of output
}

/*
 * Try to open the given bcl/scl file.
 * If that doesn't work, try appending ".gz" and gzopen it
 * If *that* doesn't work, return the error message in the errmsg field
 */
bclfile_t *bclfile_open(char *fname)
{
    bclfile_t *bclfile = calloc(1, sizeof(bclfile_t));
    bclfile->current_cluster = 0;
    bclfile->total_clusters = 0;
    bclfile->gzhandle = NULL;
    bclfile->file_type = BCL_UNKNOWN;
    bclfile->current_base = 0;
    bclfile->filename = strdup(fname);
    bclfile->nbins = 0;
    bclfile->qbin = ia_init(10);
    bclfile->qscore = ia_init(10);
    bclfile->tiles = va_init(500,free);
    bclfile->current_block = NULL;
    bclfile->current_block_size = 0;
    char *gzfname = NULL;
    int r, n;

    // need to find if this is a BCL or SCL file
    if (strstr(fname,".bcl")) bclfile->file_type = BCL_BCL;
    else if (strstr(fname,".scl")) bclfile->file_type = BCL_SCL;
    else if (strstr(fname,".cbcl")) bclfile->file_type = BCL_CBCL;

    if (bclfile->file_type == BCL_UNKNOWN) {
        fprintf(stderr,"BCL file '%s' is of unknown type\n", fname);
        bclfile_close(bclfile); return NULL;
    }

    // Actually open the file. Try fname, fname.gz, fname.bgz

    bclfile->fhandle = open(fname, O_RDONLY);
    if (bclfile->fhandle == -1) {
        gzfname = calloc(1,strlen(fname)+8);
        if (!gzfname) {
            fprintf(stderr,"Failed to allocated memory for gzfname");
            return NULL;
        }
        strcpy(gzfname,fname); strcat(gzfname,".gz");
        bclfile->gzhandle = gzopen(gzfname,"r");
        if (bclfile->gzhandle == NULL) {
            strcpy(gzfname,fname); strcat(gzfname,".bgzf");
            bclfile->gzhandle = gzopen(gzfname,"r");
            if (bclfile->gzhandle == NULL) {
                bclfile->errmsg = strdup(strerror(errno));
                free(gzfname);
                return bclfile;
            }
        }
    }

    // File is open. Read and parse header.

    while (true) {
        if (bclfile->file_type == BCL_BCL || bclfile->file_type == BCL_SCL) {
            if (bclfile->gzhandle) r = gzread(bclfile->gzhandle,(void *)&bclfile->total_clusters,4);
            else                   r = read(bclfile->fhandle, (void *)&bclfile->total_clusters, 4);
            if (r <= 0) break;
        }

        if (bclfile->file_type == BCL_CBCL) {
            r = read(bclfile->fhandle, (void *)&bclfile->version, sizeof(bclfile->version)); if (r<=0) break;
            r = read(bclfile->fhandle, (void *)&bclfile->header_size, sizeof(bclfile->header_size)); if (r<=0) break;
            r = read(bclfile->fhandle, (void *)&bclfile->bits_per_base, sizeof(bclfile->bits_per_base)); if (r<=0) break;
            r = read(bclfile->fhandle, (void *)&bclfile->bits_per_qual, sizeof(bclfile->bits_per_qual)); if (r<=0) break;
            r = read(bclfile->fhandle, (void *)&bclfile->nbins, sizeof(bclfile->nbins)); if (r<=0) break;
            for (n=0; n < bclfile->nbins; n++) {
                uint32_t qbin, qscore;
                r = read(bclfile->fhandle, (void *)&qbin, sizeof(qbin)); if (r<=0) break;
                r = read(bclfile->fhandle, (void *)&qscore, sizeof(qscore)); if (r<=0) break;
                ia_push(bclfile->qbin, qbin);
                ia_push(bclfile->qscore, qscore);
            }
            r = read(bclfile->fhandle, (void *)&bclfile->ntiles, sizeof(bclfile->ntiles)); if (r<=0) break;
            for (n=0; n < bclfile->ntiles; n++) {
                tilerec_t *tilerec = calloc(1, sizeof(tilerec_t));
                uint32_t x;
                r = read(bclfile->fhandle, (void *)&x, sizeof(x)); if (r<=0) break;
                tilerec->tilenum = x;
                r = read(bclfile->fhandle, (void *)&x, sizeof(x)); if (r<=0) break;
                tilerec->nclusters = x;
                r = read(bclfile->fhandle, (void *)&x, sizeof(x)); if (r<=0) break;
                tilerec->uncompressed_blocksize = x;
                r = read(bclfile->fhandle, (void *)&x, sizeof(x)); if (r<=0) break;
                tilerec->compressed_blocksize = x;
                va_push(bclfile->tiles, tilerec);
            }
            r = read(bclfile->fhandle, (void *)&bclfile->pfFlag, sizeof(bclfile->pfFlag));
        }
        break;
    }

    if (r <= 0) {
        fprintf(stderr,"failed to read header from bcl file '%s'\n", bclfile->gzhandle ? gzfname : fname);
        bclfile_close(bclfile); bclfile = NULL;
    }

    if (bclfile->file_type == BCL_CBCL) {
        if (bclfile->bits_per_base != 2) {
            fprintf(stderr,"CBCL file '%s' has bits_per_base %d : expecting 2\n", (bclfile->gzhandle ? gzfname : fname), bclfile->bits_per_base);
            bclfile_close(bclfile); bclfile = NULL;
        }
        if (bclfile->bits_per_qual != 2) {
            fprintf(stderr,"CBCL file '%s' has bits_per_qual %d : expecting 2\n", (bclfile->gzhandle ? gzfname : fname), bclfile->bits_per_qual);
            bclfile_close(bclfile); bclfile = NULL;
        }
    }

    free(gzfname);
    return bclfile;
}

void bclfile_seek(bclfile_t *bcl, int cluster)
{
    if (bcl->gzhandle) {
        gzseek(bcl->gzhandle, (z_off_t)(4 + cluster), SEEK_SET);
    } else {
        lseek(bcl->fhandle, (off_t)(4 + cluster), SEEK_SET);
    }
}

int bclfile_seek_tile(bclfile_t *bcl, int tile)
{
    int offset;
    bool found = false;
    tilerec_t *ti;
    char *compressed_block = NULL;
    int r;

    if (bcl->file_type != BCL_CBCL) {
        fprintf(stderr,"ERROR: calling bcl_tile_seek() for non CBCL file type\n");
        return -1;
    }

    // If the tile is not in this CBCL file, it's not an error, we just ignore the request
    if (bcl->surface != bcl_tile2surface(tile)) {
        return 0;
    }

    offset = bcl->header_size;
    for (int n=0; n < bcl->tiles->end; n++) {
        ti = (tilerec_t *)bcl->tiles->entries[n];
        if (ti->tilenum == tile) {
            found = true;
            break;
        }
        offset += ti->compressed_blocksize;
    }
    if (!found) {
        fprintf(stderr,"bclfile_seek_tile(%d) : no such tile\n", tile);
        return -1;
    }

    bcl->current_tile = ti;
    lseek(bcl->fhandle, (off_t)offset, SEEK_SET);
    free(bcl->current_block);
    bcl->current_block_size = ti->uncompressed_blocksize;
    bcl->current_block = malloc(ti->uncompressed_blocksize);
    if (!bcl->current_block) {
        fprintf(stderr,"bclfile_seek_tile(%d): failed to malloc current_block\n", tile);
        return -1;
    }
    compressed_block = malloc(ti->compressed_blocksize);
    if (!compressed_block) {
        fprintf(stderr,"bclfile_seek_tile(%d): failed to malloc compressed_block\n", tile);
        return -1;
    }
    r = read(bcl->fhandle, (void *)compressed_block, ti->compressed_blocksize);
    if (r != ti->compressed_blocksize) {
        fprintf(stderr,"bclfile_seek_tile(%d): failed to read block: returned %d\n", tile, r);
        return -1;
    }
    r=uncompressBlock(compressed_block, ti->compressed_blocksize, bcl->current_block, ti->uncompressed_blocksize);
    free(compressed_block);
    if (r<0) {
        fprintf(stderr,"uncompressBlock() somehow failed in bclfile_seek_tile(%d)\n", tile);
        fprintf(stderr,"compressed_blocksize %d   uncompressed_blocksize %d\n", ti->compressed_blocksize, ti->uncompressed_blocksize);
        fprintf(stderr,"file: %s\nsurface %d\n", bcl->filename, bcl->surface);
    }
    return r;
}

void bclfile_close(bclfile_t *bclfile)
{
    if (bclfile->gzhandle) {
        gzclose(bclfile->gzhandle);
    } else {
        if (bclfile->fhandle != -1) close(bclfile->fhandle);
    }
    free(bclfile->filename);
    free(bclfile->errmsg);
    ia_free(bclfile->qbin);
    ia_free(bclfile->qscore);
    va_free(bclfile->tiles);
    free(bclfile->current_block);
    free(bclfile);
}

int bclfile_next(bclfile_t *bcl)
{
    int i=0;
    unsigned char c = 0;

    if (bcl->file_type == BCL_CBCL) {
        if (bcl->current_block == NULL) {
            tilerec_t *t = bcl->tiles->entries[0];
            bclfile_seek_tile(bcl, t->tilenum);
            bcl->block_index = 0;
        }
    }

    if (bcl->current_base == 0) {
        if (bcl->gzhandle) {
            i = gzgetc(bcl->gzhandle);
            if (i<0) return i;
            bcl->current_byte = i;
        } else {
            if (bcl->file_type == BCL_CBCL) {
                if (bcl->block_index >= bcl->current_block_size) {
                    fprintf(stderr,"bclfile_next(): block_index: %d   current_block_size: %d\n", bcl->block_index, bcl->current_block_size);
                    return -1;
                }
                bcl->current_byte = bcl->current_block[bcl->block_index++];
            } else {
                if (read(bcl->fhandle, (void *)&(bcl->current_byte), 1) != 1) return -1;
            }
        }
    }

    c = bcl->current_byte;

    if (bcl->file_type == BCL_CBCL) {
        int baseIndex = 0;
        switch (bcl->current_base) {
            case 1: baseIndex = (c >> 4) & 0x03;    
                    bcl->quality = (c >> 6) & 0x03;
                    break;
            case 0: baseIndex = c & 0x03;           
                    bcl->quality = (c >> 2) & 0x03;
                    break;
        }
        bcl->base = BCL_BASE_ARRAY[baseIndex];
        // convert quality bin number to quality score
        for (int n=0; n < bcl->qbin->end; n++) {
            if (bcl->quality == bcl->qbin->entries[n]) {
                bcl->quality = bcl->qscore->entries[n];
                break;
            }
        }
        if (!bcl->quality) {
            bcl->base = BCL_UNKNOWN_BASE;
        }
        bcl->current_base++;
        if (bcl->current_base > 1) bcl->current_base = 0;
    }

    if (bcl->file_type == BCL_SCL) {
        int baseIndex = 0;
        switch (bcl->current_base) {
            case 0: baseIndex = (c >> 6) & 0x03;    break;
            case 1: baseIndex = (c >> 4) & 0x03;    break;
            case 2: baseIndex = (c >> 2) & 0x03;    break;
            case 3: baseIndex = c & 0x03;           break;
        }
        bcl->base = BCL_BASE_ARRAY[baseIndex];
        bcl->current_base++;
        if (bcl->current_base > 3) bcl->current_base = 0;
    }

    if (bcl->file_type == BCL_BCL) {
        int baseIndex = c & 0x03;   // last two bits
        bcl->quality = (c & 0xfc) >> 2;     // rest of bits
        if (bcl->quality) {
            bcl->base = BCL_BASE_ARRAY[baseIndex];
        } else {
            bcl->base = BCL_UNKNOWN_BASE;
        }
    }

    if (bcl->current_base == 0) bcl->current_cluster++;
    return 0;
}

