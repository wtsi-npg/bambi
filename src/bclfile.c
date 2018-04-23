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
#include <htslib/hts_endian.h>

#include "bclfile.h"

// posix_fadvise control for NovaSeq
// bit 0 set => Tell the filesystem which tile we want next
// bit 1 set => Tell the filesystem when we've finished with the tile data
#ifndef USE_POSIX_FADVISE
#define USE_POSIX_FADVISE 3
#endif

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
    zInfo.next_in= (unsigned char *) abSrc;
    zInfo.next_out= (unsigned char *) abDst;

    int nErr, nRet= -1;
    nErr= inflateInit2( &zInfo, 15+32 );               // zlib function
    if (nErr != Z_OK) fprintf(stderr,"inflateInit() failed: %d\n", nErr);
    if ( nErr == Z_OK ) {
        nErr= inflate( &zInfo, Z_FINISH );     // zlib function
        if ( (nErr == Z_STREAM_END) || (nErr == Z_BUF_ERROR) ) {
            nRet= zInfo.total_out;
            if (nRet != nLenDst) fprintf(stderr,"inflate() returned %d: expected %d\n",nRet,nLenDst);
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


bclfile_t *bclfile_init(void)
{
    bclfile_t *bclfile = calloc(1, sizeof(bclfile_t));
    bclfile->base_ptr = 0;
    bclfile->total_clusters = 0;
    bclfile->fhandle = NULL;
    bclfile->gzhandle = NULL;
    bclfile->is_cached = 0;
    bclfile->machine_type = MT_UNKNOWN;
    bclfile->current_base = 0;
    bclfile->filename = NULL;
    bclfile->nbins = 0;
    bclfile->tiles = va_init(500,free);
    bclfile->current_block = NULL;
    bclfile->current_block_size = 0;
    bclfile->surface = 1;
    bclfile->fails = 0;
    return bclfile;
}

 /*
 * Try to open the given bcl file.
 */
static void _bclfile_open_miseq(bclfile_t *bcl)
{
    int r;
    
    bcl->fhandle = fopen(bcl->filename, "rb");
    if (bcl->fhandle == NULL) die("Can't open BCL file %s\n", bcl->filename);

    r = fread(&bcl->total_clusters, 1, 4, bcl->fhandle);
    if (r != 4) die("failed to read header from bcl file '%s'\n", bcl->filename);

    char *buffer = malloc(bcl->total_clusters);
    if (!buffer) die("Can't malloc buffer %d for file %s\n", bcl->total_clusters, bcl->filename);
    r = fread(buffer, 1, bcl->total_clusters, bcl->fhandle);
    if (r != bcl->total_clusters) die("failed to read buffer from bcl file '%s'\n", bcl->filename);
    
    free(bcl->bases); bcl->bases = malloc(bcl->total_clusters);
    free(bcl->quals); bcl->quals = malloc(bcl->total_clusters);

    for (int n=0; n < bcl->total_clusters; n++) {
        char c = buffer[n];
        int baseIndex = c & 0x03;   // last two bits
        bcl->quals[n] = (c & 0xfc) >> 2;     // rest of bits
        if (bcl->quals[n]) bcl->bases[n] = BCL_BASE_ARRAY[baseIndex];
        else               bcl->bases[n] = BCL_UNKNOWN_BASE;
    }
    free(buffer);
    bcl->base_ptr = 0;
    bcl->bases_size = bcl->total_clusters;
}

static void _bclfile_open_hiseqx(bclfile_t *bcl)
{
    int r;
    
    bcl->gzhandle = gzopen(bcl->filename, "r");
    if (!bcl->gzhandle) die("Can't open BCL file %s\n", bcl->filename);

    r = gzread(bcl->gzhandle, (void *)&bcl->total_clusters, 4);
    if (r != 4) die("failed to read header from bcl file '%s'\n", bcl->filename);

    char *buffer = malloc(bcl->total_clusters);
    if (!buffer) die("Can't malloc buffer %d for file %s\n", bcl->total_clusters, bcl->filename);
    r = gzread(bcl->gzhandle, buffer, bcl->total_clusters);
    if (r != bcl->total_clusters) die("failed to read buffer from bcl file '%s'\n", bcl->filename);
    
    free(bcl->bases); bcl->bases = malloc(bcl->total_clusters);
    free(bcl->quals); bcl->quals = malloc(bcl->total_clusters);

    for (int n=0; n < bcl->total_clusters; n++) {
        char c = buffer[n];
        int baseIndex = c & 0x03;   // last two bits
        bcl->quals[n] = (c & 0xfc) >> 2;     // rest of bits
        if (bcl->quals[n]) bcl->bases[n] = BCL_BASE_ARRAY[baseIndex];
        else               bcl->bases[n] = BCL_UNKNOWN_BASE;
    }
    free(buffer);
    bcl->base_ptr = 0;
    bcl->bases_size = bcl->total_clusters;
}

static void _bclfile_open_nextseq(bclfile_t *bcl)
{
    fprintf(stderr,"WARNING: NextSeq files have not been tested properly. Trying to open %s\n", bcl->filename);
    _bclfile_open_hiseqx(bcl);
}

static off_t find_tile_offset(bclfile_t *bcl, int tile, tilerec_t **ti_out)
{
    off_t offset = bcl->header_size;
    bool found = false;

    for (int n=0; n < bcl->tiles->end; n++) {
        tilerec_t *ti = (tilerec_t *)bcl->tiles->entries[n];
        if (ti->tilenum == tile) {
            found = true;
            if (ti_out != NULL) *ti_out = ti;
            break;
        }
        offset += ti->compressed_blocksize;
    }
    return found ? offset : -1;
}

static void _bclfile_open_novaseq(bclfile_t *bclfile, int tile)
{
    int r;
    uint32_t n, m;
    uint8_t buffer[4096];
    const uint32_t qbins_in_buffer = sizeof(buffer) / (2 * 4);
    const uint32_t tiles_in_buffer = sizeof(buffer) / (4 * 4);

    bclfile->fhandle = fopen(bclfile->filename, "rb");
    if (bclfile->fhandle == NULL) {
        die("Can't open BCL file %s\n", bclfile->filename);
    }
#if (USE_POSIX_FADVISE > 0)
    // posix_fadvise() and buffered I/O don't go well together, so switch
    // to unbuffered mode.  Makes reading the header slightly less efficient
    // but shouldn't matter for the rest of the file.
    setvbuf(bclfile->fhandle, NULL, _IONBF, 0);
#endif
    // File is open. Read and parse header.

    r = fread(buffer, 2+4+1+1+4, 1, bclfile->fhandle);
    if (r != 1) goto fail;
    bclfile->version = le_to_u16(buffer);
    bclfile->header_size = le_to_u32(buffer + 2);
    bclfile->bits_per_base = buffer[6];
    bclfile->bits_per_qual = buffer[7];
    bclfile->nbins = le_to_u32(buffer + 8);
    for (n = 0; n < bclfile->nbins;) {
        uint32_t last_n = bclfile->nbins < n + qbins_in_buffer ? bclfile->nbins : n + qbins_in_buffer;
        r = fread(buffer, 8, last_n - n, bclfile->fhandle);
        if (r != last_n - n) goto fail;
        for (m = 0; n < last_n; n++, m += 8) {
            uint32_t qbin = le_to_u32(buffer + m);
            uint32_t qscore = le_to_u32(buffer + m + 4);
            bclfile->qbin[qbin] = qscore;
        }
    }
    r = fread(&bclfile->ntiles, sizeof(bclfile->ntiles), 1, bclfile->fhandle);
    if (r!=1) goto fail;
    for (n = 0; n < bclfile->ntiles; ) {
        uint32_t last_n = bclfile->ntiles < n + tiles_in_buffer ? bclfile->ntiles : n + tiles_in_buffer;
        r = fread(buffer, 16, last_n - n, bclfile->fhandle);
        if (r != last_n - n) goto fail;
        for (m = 0; n < last_n; m += 16, n++) {
            tilerec_t *tilerec = calloc(1, sizeof(tilerec_t));
            if (!tilerec) die("Out of memory\n");
            tilerec->tilenum = le_to_u32(buffer + m);
            tilerec->nclusters = le_to_u32(buffer + m + 4);
            tilerec->uncompressed_blocksize = le_to_u32(buffer + m + 8);
            tilerec->compressed_blocksize = le_to_u32(buffer + m + 12);
            va_push(bclfile->tiles, tilerec);
            if (!bclfile->current_tile) bclfile->current_tile = tilerec;
        }
    }
    r = fread(&bclfile->pfFlag, sizeof(bclfile->pfFlag), 1, bclfile->fhandle);
    if (r!=1) goto fail;

    if (bclfile->bits_per_base != 2) {
        die("CBCL file '%s' has bits_per_base %d : expecting 2\n", (bclfile->filename), bclfile->bits_per_base);
    }
    if (bclfile->bits_per_qual != 2) {
        die("CBCL file '%s' has bits_per_qual %d : expecting 2\n", (bclfile->filename), bclfile->bits_per_qual);
    }
#if (USE_POSIX_FADVISE & 1) > 0
    if (tile >= 0) {
        tilerec_t *ti = NULL;
        off_t offset = find_tile_offset(bclfile, tile, &ti);
        if (offset >= 0) {
            posix_fadvise(fileno(bclfile->fhandle), offset, ti->compressed_blocksize, POSIX_FADV_WILLNEED);
        }
    }
#endif
    return;
 fail:
    die("failed to read header from bcl file '%s'\n", bclfile->filename);
}

bclfile_t *bclfile_open(char *fname, MACHINE_TYPE mt, int tile)
{
    bclfile_t *bclfile = bclfile_init();
    bclfile->filename = strdup(fname);
    bclfile->machine_type = mt;

    switch(mt) {
        case MT_MISEQ: _bclfile_open_miseq(bclfile); break;
        case MT_NEXTSEQ: _bclfile_open_nextseq(bclfile); break;
        case MT_HISEQX: _bclfile_open_hiseqx(bclfile); break;
        case MT_NOVASEQ: _bclfile_open_novaseq(bclfile, tile); break;
        default: die("Unknown machine type\n");
    }
    return bclfile;
}



int bclfile_seek_cluster(bclfile_t *bcl, int cluster)
{
    if (bcl->gzhandle) {
        if (gzseek(bcl->gzhandle, (z_off_t)(4 + cluster), SEEK_SET) < 0) {
            int e;
            die("gzseek failed: %s\n", gzerror(bcl->gzhandle, &e));
        }
    } else {
        if (fseeko(bcl->fhandle, (off_t)(4 + cluster), SEEK_SET) < 0) {
            die("fseeko failed: %s", strerror(errno));
        }
    }
    return 0;
}

int bclfile_seek_tile(bclfile_t *bcl, int tile, filter_t *filter, int next_tile)
{
    off_t offset;
    tilerec_t *ti = NULL;
    char *compressed_block = NULL;
    char *uncompressed_block = NULL;
    int r;

    if (bcl->machine_type != MT_NOVASEQ) {
        fprintf(stderr,"ERROR: calling bcl_tile_seek() for non CBCL file type\n");
        return -1;
    }

    // If the tile is not in this CBCL file, it's not an error, we just ignore the request
    if (bcl->surface != bcl_tile2surface(tile)) {
        return 0;
    }

    // First, find the correct tile in the tile list
    offset = find_tile_offset(bcl, tile, &ti);
    if (offset < 0) {
        fprintf(stderr,"bclfile_seek_tile(%d) : no such tile\n", tile);
        return -1;
    }

    bcl->current_tile = ti;

    // Read and uncompress the record for this tile
    if (fseeko(bcl->fhandle, offset, SEEK_SET) < 0) {
        die("Couldn't seek: %s\n", strerror(errno));
    }
    uncompressed_block = malloc(ti->uncompressed_blocksize);
    if (!uncompressed_block) {
        fprintf(stderr,"bclfile_seek_tile(%d): failed to malloc uncompressed_block\n", tile);
        return -1;
    }
    compressed_block = malloc(ti->compressed_blocksize);
    if (!compressed_block) {
        fprintf(stderr,"bclfile_seek_tile(%d): failed to malloc compressed_block\n", tile);
        return -1;
    }
    r = fread(compressed_block, 1, ti->compressed_blocksize, bcl->fhandle);
    if (r != ti->compressed_blocksize) {
        fprintf(stderr,"bclfile_seek_tile(%d): failed to read block: returned %d\n", tile, r);
        return -1;
    }

#if (USE_POSIX_FADVISE & 2) > 0
    posix_fadvise(fileno(bcl->fhandle), offset, ti->compressed_blocksize, POSIX_FADV_DONTNEED);
#endif
#if (USE_POSIX_FADVISE & 1) > 0
    if (next_tile >= 0) {
        tilerec_t *next_ti = NULL;
        off_t next_offset = find_tile_offset(bcl, next_tile, &next_ti);
        if (next_offset >= 0) {
            posix_fadvise(fileno(bcl->fhandle), next_offset, next_ti->compressed_blocksize, POSIX_FADV_WILLNEED);
        }
    }
#endif

    r=uncompressBlock(compressed_block, ti->compressed_blocksize, uncompressed_block, ti->uncompressed_blocksize);
    free(compressed_block);
    if (r<0) {
        fprintf(stderr,"uncompressBlock() somehow failed in bclfile_seek_tile(%d)\n", tile);
        fprintf(stderr,"compressed_blocksize %d   uncompressed_blocksize %d\n", ti->compressed_blocksize, ti->uncompressed_blocksize);
        fprintf(stderr,"file: %s\nsurface %d\n", bcl->filename, bcl->surface);
        return r;
    }

    bcl->bases_size = ti->uncompressed_blocksize * 2;   // NovaSeq stores 2 bases and 2 quals per byte
    free(bcl->bases); bcl->bases = malloc(bcl->bases_size);
    free(bcl->quals); bcl->quals = malloc(bcl->bases_size);
    if (!bcl->bases || !bcl->quals) { fprintf(stderr,"Can't malloc memory for bases or quals in bclfile_seek_tile()"); return -1; }
    int b=0;
    for (int n=0, u=0; n < ti->uncompressed_blocksize; n++, u+=2) {
        char c = uncompressed_block[n];
        int baseIndex = 0, qbin = 0;
        unsigned char base, qscore=0;
        baseIndex = c & 0x03;           
        qbin = (c >> 2) & 0x03;
        qscore = bcl->qbin[qbin];
        if (qscore) base = BCL_BASE_ARRAY[baseIndex];
        else        base = BCL_UNKNOWN_BASE;
        if (!filter || bcl->pfFlag || (filter->buffer[u] & 0x01)) {
            bcl->bases[b] = base;
            bcl->quals[b] = qscore;
            b++;
        }

        baseIndex = (c >> 4) & 0x03;    
        qbin = (c >> 6) & 0x03;
        qscore = bcl->qbin[qbin];
        if (qscore) base = BCL_BASE_ARRAY[baseIndex];
        else        base = BCL_UNKNOWN_BASE;
        if (!filter || bcl->pfFlag || (filter->buffer[u+1] & 0x01)) {
            bcl->bases[b] = base;
            bcl->quals[b] = qscore;
            b++;
        }
    }
    bcl->base_ptr = 0;
    free(uncompressed_block);
    return 0;
}

void bclfile_close(bclfile_t *bclfile)
{
    if (bclfile->is_cached) return;
    if (bclfile->gzhandle) if (gzclose(bclfile->gzhandle) != Z_OK) die("Couldn't gzclose BCL file [%s]\n", bclfile->filename);
    if (bclfile->fhandle != NULL) if (fclose(bclfile->fhandle)) die("Couldn't close BCL file [%s] Handle %d\n", bclfile->filename, bclfile->fhandle);
    free(bclfile->filename);
    free(bclfile->errmsg);
    va_free(bclfile->tiles);
    free(bclfile->current_block);
    free(bclfile->bases);
    free(bclfile->quals);
    free(bclfile);
}

char bclfile_base(bclfile_t *bcl, int cluster)
{
    if (cluster >= bcl->bases_size) die("Cluster %d greater than %d in BCL file %s\n", cluster, bcl->bases_size, bcl->filename);
    return bcl->bases[cluster];
}

int bclfile_quality(bclfile_t *bcl, int cluster)
{
    if (cluster >= bcl->bases_size) die("Cluster %d greater than %d in BCL file %s\n", cluster, bcl->bases_size, bcl->filename);
    return bcl->quals[cluster];
}

int bclfile_load_tile(bclfile_t *bcl, int tile, filter_t *filter, int next_tile)
{
    int retval = 1;

    if (bcl->machine_type == MT_NOVASEQ) retval = bclfile_seek_tile(bcl, tile, filter, next_tile);
    else if (bcl->machine_type == MT_NEXTSEQ) retval = bclfile_seek_cluster(bcl, tile);

    return retval;
}


