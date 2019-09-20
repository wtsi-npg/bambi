/*  update.c -- Make changes to a BAM file.

    Copyright (C) 2019 Genome Research Ltd.

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
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <regex.h>
#include <inttypes.h>

#include "bamit.h"
#include "hash_table.h"

/*
 * structure to hold options
 */
typedef struct {
    char *input_name;
    char *output_name;
    bool verbose;
    char *argv_list;
    char *input_fmt;
    char *output_fmt;
} opts_t;

static void free_opts(opts_t* opts)
{
    if (!opts) return;
    free(opts->input_name);
    free(opts->output_name);
    free(opts->argv_list);
    free(opts->input_fmt);
    free(opts->output_fmt);
    free(opts);
}


/*
 * display usage information
 */
static void usage(FILE *write_to)
{
    fprintf(write_to,
"Usage: bambi update [options] <infile> <outfile>\n"
"\n"
"Options:\n"
"  -v   --verbose                       verbose output\n"
"       --input-fmt                     format of input file [sam/bam/cram]\n"
"       --output-fmt                    format of output file [sam/bam/cram]\n"
"\n"
"<infile> defaults to STDIN\n"
"<outfile> defaults to STDOUT\n"
);
}

/*
 * Takes the command line options and turns them into something we can understand
 */
static opts_t* parse_args(int argc, char *argv[])
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char* optstring = "v";

    static const struct option lopts[] = {
        { "verbose",                    0, 0, 'v' },
        { "input-fmt",                  1, 0, 0 },
        { "output-fmt",                 1, 0, 0 },
        { NULL, 0, NULL, 0 }
    };

    opts_t* opts = calloc(sizeof(opts_t), 1);
    if (!opts) { perror("cannot allocate option parsing memory"); return NULL; }

    opts->argv_list = stringify_argv(argc+1, argv-1);
    if (opts->argv_list[strlen(opts->argv_list)-1] == ' ') opts->argv_list[strlen(opts->argv_list)-1] = 0;

    // set defaults
    opts->verbose = false;

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, optstring, lopts, &option_index)) != -1) {
        const char *arg;
        switch (opt) {
        case 'v':   opts->verbose = true;
                    break;
        case 0:     arg = lopts[option_index].name;
                         if (strcmp(arg, "input-fmt") == 0)     opts->input_fmt = strdup(optarg);
                    else if (strcmp(arg, "output-fmt") == 0)    opts->output_fmt = strdup(optarg);
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
    if (argc > 1) opts->output_name = strdup(argv[1]);
    optind = 0;

    // some validation and tidying

    // input defaults to stdin
    if (!opts->input_name) opts->input_name = strdup("-");

    // output defaults to stdout
    if (!opts->output_name) opts->output_name = strdup("-");

    return opts;
}

/*
 * stolen from htslib sam.c
 */
static int do_realloc_bam_data(bam1_t *b, size_t desired)
{
    uint32_t new_m_data;
    uint8_t *new_data;
    new_m_data = desired;
    kroundup32(new_m_data);
    if (new_m_data < desired) {
        return -1;
    }
    new_data = realloc(b->data, new_m_data);
    if (!new_data) return -1;
    b->data = new_data;
    b->m_data = new_m_data;
    return 0;
}

static inline int realloc_bam_data(bam1_t *b, size_t desired)
{
    if (desired <= b->m_data) return 0;
    return do_realloc_bam_data(b, desired);
}

static inline int possibly_expand_bam_data(bam1_t *b, size_t bytes) {
    uint32_t new_len = b->l_data + bytes;

    if (new_len > INT32_MAX || new_len < b->l_data) {
        return -1;
    }
    if (new_len <= b->m_data) return 0;
    return do_realloc_bam_data(b, new_len);
}
/*
 * end of htslib sam.c
 */

//
// Set qname
//
// This *really* needs to be moved into htslib
//
static void set_qname(bam1_t *rec, char *qname)
{
    int old_len = rec->core.l_qname;
    int new_len = strlen(qname) + 1;

    int extranul = (new_len%4 != 0) ? (4 - new_len%4) : 0;

    int new_data_len = rec->l_data - old_len + new_len + extranul;
    possibly_expand_bam_data(rec, new_data_len - rec->l_data);

    uint8_t *new_data = calloc(1, new_data_len);
    memcpy(new_data, qname, new_len);
    uint8_t *p = new_data + new_len;
    for (int n=0; n<extranul; n++) *p++=0;
    memcpy(p, rec->data + rec->core.l_qname, rec->l_data - rec->core.l_qname);

    free(rec->data); rec->data=new_data;
    rec->l_data = new_data_len;
    rec->core.l_qname = new_len + extranul;
    rec->core.l_extranul = extranul;

}

static int countChars(char *str, char c)
{
    int n=0;
    while (*str) n += ((*str++) == c) ? 1 : 0;
    return n;
}

/*
 * process one BAM record, and store accumulated results in 'results'
 */
static int updateRecord(bam1_t *rec, BAMit_t *bam_out, opts_t *opts, HashTable *RGhash)
{
    char *qname = bam_get_qname(rec);
    uint8_t *tag;
    char *rgid = NULL;
    char *f1=NULL, *f2=NULL, *f3=NULL;

    // does it look like it doesn't need a prefix?
    if (countChars(qname, ':') > 3) return 0;

    // look up the RG tag
    tag = bam_aux_get(rec, "RG"); if (tag) rgid = bam_aux2Z(tag);
    if (!rgid) die("no RG tag found for record %s\n", qname);

    HashItem *hi = HashTableSearch(RGhash, rgid, strlen(rgid));
    if (!hi) die("Can't find RG tag in hash for record %s\n", qname);

    char *pu = strdup(hi->data.p);
    char *new_qname = calloc(1, strlen(qname) + strlen(pu) + 2);

    char *saveptr;
    char *s = strtok_r(pu, "_", &saveptr);    
    s = strtok_r(NULL, "_", &saveptr); if (s) f1 = strdup(s);
    s = strtok_r(NULL, "_", &saveptr); if (s) f2 = strdup(s);
    s = strtok_r(NULL, "_", &saveptr); if (s) f3 = strdup(s);

    if (!f1 || !f2 || !f3) die("PU is in incorrect format: %s\n", pu);

    strcpy(new_qname, f1);
    strcat(new_qname, ":");
    strcat(new_qname, f2);
    strcat(new_qname, ":");
    strcat(new_qname, f3);
    strcat(new_qname, ":");
    strcat(new_qname, qname);

    set_qname(rec, new_qname);

    free(pu); free(new_qname); free(f1); free(f2); free(f3);
    return 0;
}

//
// Create a hash of @RG header lines where key is RG ID 
// and value is the PU field of the line
//
static HashTable *loadRGhash(sam_hdr_t *hdr)
{
    HashTable *RGhash = HashTableCreate(0, HASH_DYNAMIC_SIZE|HASH_FUNC_JENKINS);
    kstring_t key;
    kstring_t val;
    ks_initialize(&key); ks_initialize(&val);

    int nRG = sam_hdr_count_lines(hdr, "RG");
    for (int n=0; n<nRG; n++) {
        int added;
        HashData data;
        if (sam_hdr_find_tag_pos(hdr, "RG", n, "ID", &key)) die("Can't find RG ID for entry %d\n", n);
        if (sam_hdr_find_tag_pos(hdr, "RG", n, "PU", &val)) die("Can't find RG PU for entry %d\n", n);
        data.p = strdup(ks_str(&val));
        HashTableAdd(RGhash, ks_str(&key), 0, data, &added);
        if (!added) die("loadRGhash(): failed to add %s : %s to hash\n", ks_str(&key), ks_str(&val));
    }
    ks_free(&key); ks_free(&val);
    return RGhash;
}

/*
 * Main code
 */
static int update(opts_t* opts)
{
    int retcode = 0;
    BAMit_t *bam_in = NULL;
    BAMit_t *bam_out = NULL;
    bam1_t *rec = NULL;
    HashTable *RGhash;

    // Open input BAM file
    bam_in = BAMit_open(opts->input_name, 'r', opts->input_fmt, 0, NULL);
    if (!bam_in) return 1;

    // Open output BAM file
    bam_out = BAMit_open(opts->output_name, 'w', opts->output_fmt, 0, NULL);
    if (!bam_out) return 1;

    // copy input to output header
    sam_hdr_destroy(bam_out->h); bam_out->h = sam_hdr_dup(bam_in->h);
    sam_hdr_add_pg(bam_out->h, "bambi",
                   "VN", bambi_version(),
                   "CL", opts->argv_list,
                   "DS", "update BAM file",
                   NULL, NULL);
    if (sam_hdr_write(bam_out->f, bam_out->h)) die("update(): Can't write header");

    RGhash = loadRGhash(bam_in->h);

    // Read and process each record in the input BAM
    while ( (rec = BAMit_next(bam_in)) ) {
        if (updateRecord(rec, bam_out, opts, RGhash)) { retcode = 1; break; }
        if (sam_write1(bam_out->f, bam_out->h, rec) < 0) die("Failed to write record\n");
    }

    // tidy up after us
    HashTableDestroy(RGhash, 1);
    BAMit_free(bam_in);
    BAMit_free(bam_out);

    return retcode;
}

/*
 * called from bambi to update one or more bam files
 *
 * Parse the command line arguments, then call the main update function
 *
 * returns 0 on success, 1 if there was a problem
 */
int main_update(int argc, char *argv[])
{
    int ret = 1;

    opts_t* opts = parse_args(argc, argv);
    if (opts) {
        ret = update(opts);
    }
    free_opts(opts);
    return ret;
}
