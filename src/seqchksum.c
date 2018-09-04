/*  seqchksum.c -- Calculate chksums for a BAM file.

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

#include <bambi.h>
#include <assert.h>
#include <htslib/sam.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <unistd.h>
#include <regex.h>
#include <htslib/khash.h>
#include <cram/sam_header.h>
#include <inttypes.h>

#include "seqchksum.h"
#include "bamit.h"
#include "hash_table.h"
#include "parse_bam.h"
#include "crc.h"
#include "hash_table.h"

#define xstr(s) str(s)
#define str(s) #s

//static char *hash_names[] = {"crc32", "crc32prod"};

HASH_TYPE decode_hash_name(char *name)
{
    HASH_TYPE hash = HASH_UNKNOWN;
    if (strcmp(name, "crc32") == 0)     hash = HASH_CRC32;
    if (strcmp(name, "crc32prod") == 0)  hash = HASH_CRC32PROD;
    return hash;
}

/*
 * structure to hold options
 */
typedef struct {
    char *input_name;
    bool verbose;
    char *argv_list;
    char *input_fmt;
    HASH_TYPE hash;
} opts_t;

static void free_opts(opts_t* opts)
{
    if (!opts) return;
    free(opts->input_name);
    free(opts->argv_list);
    free(opts->input_fmt);
    free(opts);
}


/*
 * display usage information
 */
static void usage(FILE *write_to)
{
    fprintf(write_to,
"Usage: bambi seqchksum [options] <filename>\n"
"\n"
"Options:\n"
"  -v   --verbose                       verbose output\n"
"       --input-fmt                     format of input file [sam/bam/cram]\n"
"       --hash                          Hash type [default: crcprod]\n"
);
}

/*
 * Takes the command line options and turns them into something we can understand
 */
static opts_t* parse_args(int argc, char *argv[])
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char* optstring = "i:vb:";

    static const struct option lopts[] = {
        { "hash",                       1, 0, 0 },
        { "input",                      1, 0, 'i' },
        { "verbose",                    0, 0, 'v' },
        { "input-fmt",                  1, 0, 0 },
        { NULL, 0, NULL, 0 }
    };

    opts_t* opts = calloc(sizeof(opts_t), 1);
    if (!opts) { perror("cannot allocate option parsing memory"); return NULL; }

    opts->argv_list = stringify_argv(argc+1, argv-1);
    if (opts->argv_list[strlen(opts->argv_list)-1] == ' ') opts->argv_list[strlen(opts->argv_list)-1] = 0;

    // set defaults
    opts->verbose = false;
    opts->hash = DEFAULT_HASH_TYPE;

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, optstring, lopts, &option_index)) != -1) {
        const char *arg;
        switch (opt) {
        case 'i':   opts->input_name = strdup(optarg);
                    break;
        case 'v':   opts->verbose = true;
                    break;
        case 0:     arg = lopts[option_index].name;
                         if (strcmp(arg, "hash") == 0)               opts->hash = decode_hash_name(optarg);
                    else if (strcmp(arg, "input-fmt") == 0)                  opts->input_fmt = strdup(optarg);
                    else {
                        printf("\nUnknown option: %s\n\n", arg); 
                        usage(stdout); free_opts(opts);
                        return NULL;
                    }
                    break;
        default:    printf("Unknown option: '%c'\n", opt);
            /* else fall-through */
        case '?':   usage(stdout); free_opts(opts); return NULL;
        }
    }

    argc -= optind;
    argv += optind;

    if (argc > 0) opts->input_name = strdup(argv[0]);
    optind = 0;

    // some validation and tidying

    // input defaults to stdin
    if (!opts->input_name) opts->input_name = strdup("-");

    // TODO: list valid hash types here
    if (opts->hash == HASH_UNKNOWN) {
        fprintf(stderr, "Unknown has type\n");
        return NULL;
    }

    return opts;
}

static void update_crc(uint32_t *crc, uint8_t *data, size_t length)
{
    *crc = crc32(*crc, data, length);
}

uint64_t update_digest_crc32prod(uint64_t digest, uint32_t val)
{
    static uint64_t const MERSENNE31 = 0x7FFFFFFFull; // Mersenne Prime 2^31 - 1
    val &= MERSENNE31;
    if (!val || val==MERSENNE31) val=1;
    digest = (digest * val)%MERSENNE31;
    return digest;
}

uint64_t update_digest_crc32(uint64_t digest, uint32_t val)
{
    return digest + val;
}

uint64_t update_digest(HASH_TYPE hash, uint64_t digest, uint32_t val)
{
    switch (hash) {
        case HASH_CRC32:        digest = update_digest_crc32(digest,val);
                                break;
        case HASH_CRC32PROD:    digest = update_digest_crc32prod(digest,val);
                                break;
        case HASH_UNKNOWN:      break;
    }
    return digest;
}

static void init_digest_line(HASH_TYPE hash, digest_line_t *d)
{
    int v = 0;
    if (hash == HASH_CRC32) v=0;
    if (hash == HASH_CRC32PROD) v=1;

    d->count[0] = 0;    d->count[1] = 0;
    for (int c=0; c<4; c++) {
        for (int p=0; p<2; p++) {
            d->chksum[c][p] = v;
        }
    }
}

static void update_digest_line(HASH_TYPE hash, bool pass, digest_line_t *dline, uint32_t crc, int c)
{
    if (!c) dline->count[0]++;
    dline->chksum[c][0] = update_digest(hash, dline->chksum[c][0], crc);
    if (pass) {
        if (!c) dline->count[1]++;
        dline->chksum[c][1] = update_digest(hash, dline->chksum[c][1], crc);
    }
}

/*
 * process one BAM record, and store accumulated results in 'results'
 */
int seqchksum_processRecord(bam1_t *rec, HASH_TYPE hash, chksum_results_t *results)
{
    uint32_t crc = 0;

    uint16_t aflags = rec->core.flag;
    uint8_t *seq = get_read(rec);
    uint8_t *qual = get_quality(rec);
    uint16_t flag_mask = BAM_FPAIRED | BAM_FREAD1 | BAM_FREAD2;
    uint8_t flags = (aflags & flag_mask) & 0xFF;
    bool pass = !(aflags & BAM_FQCFAIL);;
    char *qname = bam_get_qname(rec);
    uint8_t *tag;
    char *rgid;
    HashItem *hi;
    HashData hd;
    int newitem;
    digest_line_t *dline_all;
    digest_line_t *dline;

    // look up the RG tag
    tag = bam_aux_get(rec, "RG");
    //hd.p = malloc(sizeof(digest_line_t));
    if (tag) rgid = bam_aux2Z(tag);
    else     rgid = "";

    hd.p = NULL;
    hi = HashTableAdd(results->rgHash, rgid, 0, hd, &newitem);
    if (newitem) {
        hi->data.p = malloc(sizeof(digest_line_t));
        dline = hi->data.p;
        init_digest_line(hash,dline);
    } else {
        dline = hi->data.p;
    }

    dline_all = &(results->all);

    // flags + sequence chksum
    update_crc(&crc,&flags,1);
    update_crc(&crc,seq,strlen((char*)seq));

    update_digest_line(hash, pass, dline, crc, 0);
    update_digest_line(hash, pass, dline_all, crc, 0);

    // flags + sequence + quality chksum (don't reset crc, just add quality)
    update_crc(&crc,qual,strlen((char*)qual));
    update_digest_line(hash, pass, dline, crc, 2);
    update_digest_line(hash, pass, dline_all, crc, 2);

    // name + flags + sequence chksum
    crc = 0;
    update_crc(&crc, (uint8_t *)qname, strlen(qname)+1);
    update_crc(&crc, &flags, 1);
    update_crc(&crc,seq,strlen((char*)seq));
    update_digest_line(hash, pass, dline, crc, 1);
    update_digest_line(hash, pass, dline_all, crc, 1);

    // flags + sequence + tags chksum
    crc = 0;
    update_crc(&crc, &flags, 1);
    update_crc(&crc,seq,strlen((char*)seq));
    tag = bam_aux_get(rec,"BC"); if (tag) update_crc(&crc,tag-2,aux_type2size(tag)+3);
    tag = bam_aux_get(rec,"FI"); if (tag) update_crc(&crc,tag-2,aux_type2size(tag)+3);
    tag = bam_aux_get(rec,"QT"); if (tag) update_crc(&crc,tag-2,aux_type2size(tag)+3);
    tag = bam_aux_get(rec,"RT"); if (tag) update_crc(&crc,tag-2,aux_type2size(tag)+3);
    tag = bam_aux_get(rec,"TC"); if (tag) update_crc(&crc,tag-2,aux_type2size(tag)+3);
    update_digest_line(hash, pass, dline, crc, 3);
    update_digest_line(hash, pass, dline_all, crc, 3);

    free(seq); free(qual);
    return 0;
}

static char *dformat(uint32_t v)
{
    static char buffer[64];
    sprintf(buffer, "0x%" PRIx32 , v);
    return buffer+2;
}

static void hputi(int n, hFILE *f)
{
    char b[64];
    sprintf(b,"%d",n);
    hputs(b,f);
}

static void print_dline(hFILE *f, char *key, digest_line_t *dline, int p)
{
    hputs(key, f); hputc('\t',f); 
    hputs((p ? "pass": "all"), f); hputc('\t',f);
    hputi(dline->count[p], f); hputc('\t',f);
    hputc('\t',f);
    hputs(dformat(dline->chksum[0][p]), f); hputc('\t',f);
    hputs(dformat(dline->chksum[1][p]), f); hputc('\t',f);
    hputs(dformat(dline->chksum[2][p]), f); hputc('\t',f);
    hputs(dformat(dline->chksum[3][p]), f); hputc('\n',f);
}

void chksum_print_results(hFILE *f, chksum_results_t *results)
{
    digest_line_t *dline = &(results->all);

    hputs("###\tset\tcount\t\tb_seq\tname_b_seq\tb_seq_qual\tb_seq_tags(BC,FI,QT,RT,TC)\n", f);

    print_dline(f, "all", dline, 0);
    print_dline(f, "all", dline, 1);

    HashIter *iter = HashTableIterCreate();
    HashItem *hi;
    while ( (hi = HashTableIterNext(results->rgHash, iter)) != NULL) {
        print_dline(f, hi->key, hi->data.p, 0);
        print_dline(f, hi->key, hi->data.p, 1);
    }
    HashTableIterDestroy(iter);
}

/*
 * Initialise results structure
 */
chksum_results_t *chksum_init_results(HASH_TYPE hash)
{
    chksum_results_t *results = malloc(sizeof(chksum_results_t));
    init_digest_line(hash, &(results->all));
    results->rgHash = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS);
    return results;
}

/*
 * free results structure
 */
void chksum_free_results(chksum_results_t *results)
{
    HashTableDestroy(results->rgHash, 1);
    free(results);
}
 
/*
 * Main code
 */
static int seqchksum(opts_t* opts)
{
    int retcode = 0;
    BAMit_t *bam_in = NULL;
    bam1_t *rec = NULL;

    /*
     * Open input BAM file
     */
    bam_in = BAMit_open(opts->input_name, 'r', opts->input_fmt, 0, NULL);
    if (!bam_in) return 1;

    // Initialise results structure
    chksum_results_t *results = chksum_init_results(opts->hash);

    // Read and process each record in the input BAM
    while ( (rec = BAMit_next(bam_in)) ) {
        // ignore secondary and supplementary records
        if (!(rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) {
            if (seqchksum_processRecord(rec, opts->hash, results)) { retcode = 1; break; }
        }
    }

    hFILE *f = hdopen(fileno(stdout),"w");
    if (!f) die("Can't open stdout");
    chksum_print_results(f, results);
    if (hclose(f)) die("Can't close stdout");

    // tidy up after us
    BAMit_free(bam_in);
    chksum_free_results(results);

    return retcode;
}

/*
 * called from bambi to perform index decoding
 *
 * Parse the command line arguments, then call the main decode function
 *
 * returns 0 on success, 1 if there was a problem
 */
int main_seqchksum(int argc, char *argv[])
{
    int ret = 1;

    opts_t* opts = parse_args(argc, argv);
    if (opts) {
        ret = seqchksum(opts);
    }
    free_opts(opts);
    return ret;
}
