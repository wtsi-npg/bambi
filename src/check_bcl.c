/* chack_bcl.c

   Program to sanity check one or more BCL files

    Copyright (C) 2018 Genome Research Ltd.

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

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <pthread.h>
#include <htslib/thread_pool.h>

#include "bclfile.h"

#define NTHREADS 16

int nFailed = 0, nPassed = 0;
pthread_mutex_t nFailed_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t nPassed_lock = PTHREAD_MUTEX_INITIALIZER;
int check_opts_verbose = 0;

static char *reverse_str(char *seq)
{
    char *t = seq, *s = seq + strlen(seq) - 1;
    char c;

    while (t < s) {
        c = *s;
        *s = *t;
        *t = c;
        t++;
        s--;
    }

    return seq;
}

static bool rev_strstr(char *str, char *ext)
{
    char *rev_str = strdup(str); rev_str = reverse_str(rev_str);
    char *rev_ext = strdup(ext); rev_ext = reverse_str(rev_ext);
    int r = strncmp(rev_str, rev_ext, strlen(rev_ext));
    free(rev_str); free(rev_ext);
    return r==0;
}

static void Usage(FILE *f)
{
    fprintf(f, "\n");
    fprintf(f, "check_bcl will perform a sanity check on a BCL (or CBCL) file, or on all such\n");
    fprintf(f, "files in a given directory.\n");
    fprintf(f, "If the -v (verbose) flag is given, then each file will be listed as it is checked,\n");
    fprintf(f, "else only files which fail the check will be listed\n\n");
    fprintf(f, "Usage:\n");
    fprintf(f, "check_bcl [-v] <directory>\n");
    fprintf(f, "or\n");
    fprintf(f, "check_bcl <bcl_file>\n");
    exit(1);
}

static int isDirectory(char *fname)
{
    struct stat buf;
    int is_directory = 0;
    if (stat(fname, &buf) == 0) {
        is_directory = S_ISDIR(buf.st_mode);
    }
    return is_directory;
}

// zlib version
static char *uncompressBlock(char* abSrc, int nLenSrc, char* abDst, int nLenDst )
{
    char *msg = NULL;
    z_stream zInfo ={0};
    zInfo.total_in=  zInfo.avail_in=  nLenSrc;
    zInfo.total_out= zInfo.avail_out= nLenDst;
    zInfo.next_in= (unsigned char *) abSrc;
    zInfo.next_out= (unsigned char *) abDst;

    if (nLenDst == 0) return 0;

    int nErr, nRet= -1;
    nErr= inflateInit2( &zInfo, 15+32 );               // zlib function
    if (nErr != Z_OK) store_msg(&msg,"inflateInit() failed: %d", nErr);
    if ( nErr == Z_OK ) {
        nErr= inflate( &zInfo, Z_FINISH );     // zlib function
        if (nErr == Z_STREAM_END) {
            nRet= zInfo.total_out;
            if (nRet != nLenDst) {
                store_msg(&msg, "inflate() returned %d: expected %d",nRet,nLenDst);
            }
        } else {
            char buff[64];
            sprintf(buff, "%02x %02x %02x %02x %02x", (uint8_t)abSrc[0], (uint8_t)abSrc[1], (uint8_t)abSrc[2], (uint8_t)abSrc[3], (uint8_t)abSrc[4]);
            store_msg(&msg, "inflate() returned %d for data %s", nErr, buff);
        }
    }
    inflateEnd( &zInfo );   // zlib function
    return msg;
}

#if 0
// libdeflate version
#include <libdeflate.h>
static char *uncompressBlock(char* abSrc, int nLenSrc, char* abDst, int nLenDst )
{
    char *msg = NULL;
    size_t actualOut;

    if (nLenDst == 0) return 0;

    int nErr;
    struct libdeflate_decompressor *compressor = libdeflate_alloc_decompressor();
    //LIBDEFLATEAPI *compressor = libdeflate_alloc_decompressor();

        nErr = libdeflate_gzip_decompress(compressor, abSrc, nLenSrc, abDst, nLenDst, &actualOut);
        if (nErr != LIBDEFLATE_SUCCESS) {
            char buff[64];
            sprintf(buff, "%02x %02x %02x %02x %02x", (uint8_t)abSrc[0], (uint8_t)abSrc[1], (uint8_t)abSrc[2], (uint8_t)abSrc[3], (uint8_t)abSrc[4]);
            store_msg(&msg, "inflate() returned %d for data %s", nErr, buff);
        }

    libdeflate_free_decompressor(compressor);
    return msg;
}
#endif

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

static char *checkTile(bclfile_t *bcl, int tilenum)
{
    char *msg = NULL;
    off_t offset;
    tilerec_t *ti = NULL;
    char *compressed_block = NULL;
    char *uncompressed_block = NULL;
    int r;

    offset = find_tile_offset(bcl, tilenum, &ti);
    if (offset < 0) {
        store_msg(&msg, "Can't find tile %d", tilenum);
        return msg;
    }
    if (fseeko(bcl->fhandle, offset, SEEK_SET) < 0) {
        store_msg(&msg, "Couldn't seek to %d", offset);
        return msg;
    }

    uncompressed_block = smalloc(ti->uncompressed_blocksize);
    compressed_block = smalloc(ti->compressed_blocksize);
    r = fread(compressed_block, 1, ti->compressed_blocksize, bcl->fhandle);
    if (r != ti->compressed_blocksize) {
        store_msg(&msg, "Can't read block: returned %d", r);
    } else {
        msg = uncompressBlock(compressed_block, ti->compressed_blocksize, uncompressed_block, ti->uncompressed_blocksize);
    }

    free(compressed_block);
    free(uncompressed_block);

    return msg;
}

static int checkBclFile(char *fname, int verbose)
{
    int tile = 0;
    int ret = 0;
    MACHINE_TYPE mt = MT_UNKNOWN;

         if (strstr(fname, ".bcl.gz")) mt = MT_HISEQX;
    else if (strstr(fname, ".cbcl")) mt = MT_NOVASEQ;
    else if (strstr(fname, ".bcl.bgzf")) mt = MT_NEXTSEQ;
    else if (strstr(fname, ".bcl")) mt = MT_MISEQ;

    bclfile_t *bcl;
    bcl = bclfile_open(fname, mt, tile);
    if (bcl->errmsg) {
        display("File: %s\t%s\n", fname, bcl->errmsg);
        ret = 1;
    }

    if (!bcl->errmsg) {
        if (verbose) {
            display("Filename     : %s\n", fname);
            display("Clusters     : %d\n", bcl->total_clusters);
            if (mt == MT_NOVASEQ) {
                display("Version      : %d\n", bcl->version);
                display("Header Size  : %d\n", bcl->header_size);
                display("Bits-per-base: %d\n", bcl->bits_per_base);
                display("Bits-per-qual: %d\n", bcl->bits_per_qual);
                display("nBins        : %d\n", bcl->nbins);
                display("Bins         : ");
                for (int n=0; n<bcl->nbins; n++) display("%d ", bcl->qbin[n]);
                display("\n");
                display("nTiles       : %d\n", bcl->ntiles);
                display("Tiles        :\n");
            }
        }
        if (mt == MT_NOVASEQ) {
            for (int n=0; n < bcl->tiles->end; n++) {
                tilerec_t *tile = bcl->tiles->entries[n];
                if (verbose) display("  %3d %6d %d\t%d\t%d\t", n, tile->tilenum, tile->nclusters, tile->uncompressed_blocksize, tile->compressed_blocksize);
                char *msg = checkTile(bcl,tile->tilenum);
                if (msg) ret = 1;
                if (verbose) {
                    if (msg) display("***FAIL***  %s\n", msg);
                    else     display("Ok\n");
                }
                free(msg);
            }
        }
    }

    if (!verbose && ret != 0) {
        display("Failed %s check: %s\n", (ret==1 ? "Header" : "Tile"), fname);
    }

    bclfile_close(bcl);
    return ret;
}

static char *bclname(char *d, char *f)
{
    char *bname = smalloc(strlen(d) + strlen(f) + 2);
    sprintf(bname, "%s/%s", d, f);
    return bname;
}

void *checkBclFileThread(void *arg)
{
    char *fullname = (char *)arg;

    int r = checkBclFile(fullname,0);
    if (r) {
        pthread_mutex_lock(&nFailed_lock);
        nFailed++;
        pthread_mutex_unlock(&nFailed_lock);
        if (check_opts_verbose) checkBclFile(fullname,1);
    } else {
        pthread_mutex_lock(&nPassed_lock);
        nPassed++;
        pthread_mutex_unlock(&nPassed_lock);
    }
    free(fullname);
    return NULL;
}

static void recurseThroughDirectory(char *dirname, hts_tpool *p, hts_tpool_process *q)
{
    struct dirent *dir;
    DIR *d = opendir(dirname);
    if (!d) die("Can't open directory: %s\n", dirname);

    while ( (dir = readdir(d)) != NULL) {
        // Look for a lane directory
        if (dir->d_name[0]=='.') continue;
        char *fullname = bclname(dirname,dir->d_name);
        if (isDirectory(fullname)) {
            recurseThroughDirectory(fullname, p, q);
            free(fullname);
        } else {
            // Look for bcl or cbcl files
            if (rev_strstr(dir->d_name, ".cbcl")   || 
                rev_strstr(dir->d_name, ".bcl")    ||
                rev_strstr(dir->d_name, ".bcl.gz") ||
                rev_strstr(dir->d_name, ".bcl.bgzf")
            ) {
                if (check_opts_verbose) display("%s\n",fullname);
                hts_tpool_dispatch(p, q, checkBclFileThread, fullname);
            }
        }
    }
    closedir(d);
}

static int checkRunFolder(char *dirname)
{

    hts_tpool *p = hts_tpool_init(NTHREADS);
    hts_tpool_process *q = hts_tpool_process_init(p, NTHREADS*2, 1);

    recurseThroughDirectory(dirname, p, q);

    hts_tpool_process_flush(q);
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    if (nFailed) display("Failed %d out of %d files\n", nFailed, nFailed+nPassed);
    return nFailed;
}


int main(int argc, char *argv[])
{
    int opt;
    int r = 0;
    if (argc < 2) Usage(stdout);

    const char* optstring = "v";

    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch (opt) {
        case 'v':   check_opts_verbose++; break;
        case 'h':
        case '?':   Usage(stdout); break;
        default:    Usage(stderr); break;
        }
    }

    if (optind >= argc) Usage(stderr);

    if (isDirectory(argv[optind])) r = checkRunFolder(argv[optind]);
    else                           r = checkBclFile(argv[optind], 1);

    return r;
}

