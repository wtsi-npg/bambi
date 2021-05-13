/*  -*- mode: c; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 8; -*- */

/*
 * Copyright (c) 2021, Genome Research Ltd (GRL).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *
 *     * Neither the name of the Genome Research Limited nor the
 *       names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior
 *       written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GRL ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL GRL BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author: Steven Leonard, Jan 2009
 *
 * This code generates a substitution analysis table
 *
 * Incorporated into Bambi: Jennifer Liddle <js10@sanger.ac.uk> May 2021
 *
 */

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <getopt.h>

#include "parse_bam.h"
#include "bambi_utils.h"
#include "bamit.h"

#define N_READS 3

#define LEN_SUBST       2
#define NUM_SUBST       16    // 4 ^ LEN_SUBST
#define LEN_CNTXT       4
#define NUM_CNTXT       256   // 4 ^ LEN_CNTXT

#define NBINS 51
#define ST_HILO_QUALITY  29.5

typedef struct {
    int         read;
    int         cycle;
    int         nbins;
    float       predictor_hilo;
    float       *predictor;
    long        *num_bases;
    long        *num_errors;
    long        *subst[NUM_SUBST];
    long        substH[NUM_SUBST];
    long        substL[NUM_SUBST];
    long        cntxtH[NUM_CNTXT];
    long        cntxtL[NUM_CNTXT];
    long        total_bases;
    long        total_errors;
    float       quality;
} SurvTable;

typedef struct {
    char *reportName;
    char *in_bam_file;
    int read_length[N_READS];
    int quiet;
    int verbose;
    char compression_level;
    char *input_fmt;
    char *output_fmt;

} opts_t;

static void free_opts(opts_t *opts)
{
    free(opts->reportName);
    free(opts->in_bam_file);
    free(opts->input_fmt);
    free(opts);
}


#define MAXNH       7
static int* lookup = NULL;

static void init_lookup(void) {
    int i;

    lookup = (int *) calloc(256, sizeof(int));
    for (i = 0; i < 256; i++)
    lookup[i] = -1;
    lookup['A'] = 0;
    lookup['a'] = 0;
    lookup['C'] = 1;
    lookup['c'] = 1;
    lookup['G'] = 2;
    lookup['g'] = 2;
    lookup['T'] = 3;
    lookup['t'] = 3;
}

int str2word(char *seq, int NH) {
    int i, word = -1;

    if (!lookup) {
        init_lookup();
    }

    if (NH > MAXNH)
        return word;

    word = 0;
    for (i = 0; i < NH; i++) {
        word <<= 2;
        word |= lookup[(int)seq[i]];
    }

    return word;
}

char *word2str(int word, int NH) {
    static char str[MAXNH+1];
    int i;

    if (!lookup) {
        init_lookup();
    }

    for (i = 0; i < NH; i++)
        str[i] = "ACGT"[(word >> (2*(NH-1)-2*i)) & 3];
    str[NH] = 0;

    return str;
}

static void initialiseSurvTable(SurvTable *st, int read, int cycle)
{
    int i, j;

    st->read  = read;
    st->cycle = cycle;

    st->nbins = NBINS;

    st->predictor_hilo = ST_HILO_QUALITY;

    st->predictor = (float *)smalloc(st->nbins * sizeof(float));
    for (i=0;i<st->nbins;i++)
        st->predictor[i] = i;

    st->num_bases  = (long *)smalloc(st->nbins * sizeof(long));
    st->num_errors = (long *)smalloc(st->nbins * sizeof(long));
    for (i=0;i<st->nbins;i++) {
        st->num_bases[i] = 0;
        st->num_errors[i] = 0;
    }

    for (j=0;j<NUM_SUBST;j++) {
        st->subst[j] = (long *)smalloc(st->nbins * sizeof(long));
        for (i=0;i<st->nbins;i++) 
            st->subst[j][i] = 0;
        st->substH[j]=0;
        st->substL[j]=0;
    }
    for (j=0;j<NUM_CNTXT;j++) {
        st->cntxtH[j]=0;
   	    st->cntxtL[j]=0;
    }

    st->total_bases = 0;
    st->total_errors = 0;

}

static void freeSurvEntry(SurvTable *st)
{
    if (st && st->nbins) {
        free(st->predictor);
        free(st->num_bases);
        free(st->num_errors);
        for (int i=0;i<NUM_SUBST;i++)
            free(st->subst[i]);
        st->nbins = 0;
    }
}

static void freeSurvTable(opts_t *opts, SurvTable **sts)
{
    int read, cycle;
    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[read]) continue;
        for(cycle=0;cycle<opts->read_length[read];cycle++)
        {
            SurvTable *st = sts[read] + cycle;
            freeSurvEntry(st);
        }
        free(sts[read]);
    }
    free(sts);
}

static void completeSurvTable(opts_t *opts, SurvTable **sts)
{
    float ssc = 1.0;

    for (int read=0; read < N_READS; read++) {
        if (NULL == sts[read]) continue;
        for (int cycle=0; cycle < opts->read_length[read]; cycle++) {
            SurvTable *st = sts[read] + cycle;
            long quality_bases = 0;
            long quality_errors = 0;

            for (int i=0; i < st->nbins; i++) {
                st->total_bases += st->num_bases[i];
                st->total_errors += st->num_errors[i];

                // bases in the first bin are called as N and explicitly get a quality of 0
                if( i == 0 )
                    continue;

                quality_bases  += st->num_bases[i];
                quality_errors += st->num_errors[i];
            }

            st->quality = -10.0 * log10((quality_errors + ssc)/(quality_bases + ssc));
        }
    }
}

static void writeReport(opts_t *opts, SurvTable **sts)
{
    FILE *fp = NULL;
    SurvTable **read_sts;
    int read, cycle, i, j;
    char p;

    // open report file
    if (opts->reportName) fp = fopen(opts->reportName, "w");
    else                  fp = stdout;

    if (!fp) die("ERROR: can't open report file %s: %s\n", opts->reportName, strerror(errno));

    /* generate read summary tables by summing over cycles */

    read_sts = (SurvTable **)smalloc(N_READS * sizeof(SurvTable *));

    for (read=0; read < N_READS; read++) {
        int read_length = opts->read_length[read];
        read_sts[read] = NULL;
        if (0 == read_length) continue;

        read_sts[read] = (SurvTable *)smalloc(sizeof(SurvTable));
        SurvTable *read_st = read_sts[read];

        initialiseSurvTable(read_st, read, -1);

        for (cycle=0; cycle < read_length; cycle++) {
            SurvTable *st = sts[read] + cycle;
            for(j=0;j<NUM_SUBST;j++)
                for(i=0;i<st->nbins;i++)
                    read_st->subst[j][i] += st->subst[j][i];
            for(j=0;j<NUM_SUBST;j++)
                read_st->substH[j] += st->substH[j];
            for(j=0;j<NUM_SUBST;j++)
                read_st->substL[j] += st->substL[j];
            for(j=0;j<NUM_CNTXT;j++)
                read_st->cntxtH[j] += st->cntxtH[j];
            for(j=0;j<NUM_CNTXT;j++)
                read_st->cntxtL[j] += st->cntxtL[j];
        }
    }

    /* substitution RC */

    fprintf(fp, "# Substitution error table. Use `grep ^SET | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and quality value, columns read, quality value followed by substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        SurvTable *st = read_sts[read];
        if( NULL == st) continue;

        for(i=0;i<st->nbins;i++)
        {
            fprintf(fp, "SET\t%d\t%.2f", read, st->predictor[i]);
            for (j=0;j<NUM_SUBST;j++)
            {
                char *subst;
                subst=word2str(j,LEN_SUBST);
                if (subst[0]==subst[1]) continue;
                fprintf(fp, "\t%s\t%ld", subst, st->subst[j][i]);
            }
            fprintf(fp, "\n");
        }
    }
    
    fprintf(fp, "# Mismatch substitutions high quality. Use `grep ^RCH | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and cycle, columns read, cycle then substitution and count for 12 substitutions\n");
    fprintf(fp, "# Followed by a single row with a total over all cycles for each read, columns are read, -1 then substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[read]) continue;
        /* cycle by cycle */
        for(cycle=0;cycle<opts->read_length[read];cycle++)
        {
            SurvTable *st = sts[read] + cycle;
            if( 0 == st->total_bases ) continue;
            fprintf(fp, "RCH\t%d\t%d", read, cycle);
            for (j=0;j<NUM_SUBST;j++)
            {
                char *subst;
                subst=word2str(j,LEN_SUBST);
                if (subst[0]==subst[1]) continue;
                fprintf(fp, "\t%s\t%ld", subst, st->substH[j]);
            }
            fprintf(fp, "\n");
        }
        /* and read summary (cycle = -1) */
        SurvTable *st = read_sts[read];
        fprintf(fp, "RCH\t%d\t%d", read, -1);
        for (j=0;j<NUM_SUBST;j++)
        {
            char *subst;
            subst=word2str(j,LEN_SUBST);
            if (subst[0]==subst[1]) continue;
            fprintf(fp, "\t%s\t%ld", subst, st->substH[j]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "# Mismatch substitutions low quality. Use `grep ^RCL | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and cycle, columns read, cycle then substitution and count for 12 substitutions\n");
    fprintf(fp, "# Followed by a single row with a total over all cycles for each read, columns are read, -1 then substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[read]) continue;
        /* cycle by cycle */
        for(cycle=0;cycle<opts->read_length[read];cycle++)
        {
            SurvTable *st = sts[read] + cycle;
            if( 0 == st->total_bases ) continue;
            fprintf(fp, "RCL\t%d\t%d", read, cycle);
            for (j=0;j<NUM_SUBST;j++)
            {
                char *subst;
                subst=word2str(j,LEN_SUBST);
                if (subst[0]==subst[1]) continue;
                fprintf(fp, "\t%s\t%ld", subst, st->substL[j]);
            }
            fprintf(fp, "\n");
        }
        /* and read summary (cycle = -1) */
        SurvTable *st = read_sts[read];
        fprintf(fp, "RCL\t%d\t%d", read, -1);
        for (j=0;j<NUM_SUBST;j++)
        {
            char *subst;
            subst=word2str(j,LEN_SUBST);
            if (subst[0]==subst[1]) continue;
            fprintf(fp, "\t%s\t%ld", subst, st->substL[j]);
        }
        fprintf(fp, "\n");
    }
    
    /* previous base PRC */

    fprintf(fp, "# Effect of previous base high quality. Use `grep ^PRCH | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and previous base, columns read then previous base+substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        /* read summary */
        SurvTable *st = read_sts[read];
        if( NULL == st) continue;
        int cntxt_off = 0;
        int len_cntxt = 3;
        int num_cntxt = 64;
        long *count = (long *)smalloc(num_cntxt * sizeof(long));
        for (j=0;j<num_cntxt;j++)
            count[j]=0;
        for (j=0;j<NUM_CNTXT;j++)
        {
            char *cntxt;
            int word;
            cntxt=word2str(j,LEN_CNTXT);
            cntxt+=cntxt_off;
            word=str2word(cntxt,len_cntxt);
            count[word] += st->cntxtH[j];
        }
        for (j=0,p=0;j<num_cntxt;j++)
        {
            char *cntxt;
            cntxt=word2str(j,len_cntxt);
            if (p != cntxt[0]) {
                p = cntxt[0];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "PRCH\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, count[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }
    
    fprintf(fp, "# Effect of previous base low quality. Use `grep ^PRCL | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and previous base, columns read then previous base+substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        /* read summary */
        SurvTable *st = read_sts[read];
        if( NULL == st) continue;
        int cntxt_off = 0;
        int len_cntxt = 3;
        int num_cntxt = 64;
        long *count = (long *)smalloc(num_cntxt * sizeof(long));
        for (j=0;j<num_cntxt;j++)
            count[j]=0;
        for (j=0;j<NUM_CNTXT;j++)
        {
            char *cntxt;
            int word;
            cntxt=word2str(j,LEN_CNTXT);
            cntxt+=cntxt_off;
            word=str2word(cntxt,len_cntxt);
            count[word] += st->cntxtL[j];
        }
        for (j=0,p=0;j<num_cntxt;j++)
        {
            char *cntxt;
            cntxt=word2str(j,len_cntxt);
            if (p != cntxt[0]) {
                p = cntxt[0];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "PRCL\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, count[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }
    
    /* previous base + next base PRCN */

    fprintf(fp, "# Effect of previous base and next base high quality. Use `grep ^PRCNH | cut -f 2-` to extract this part\n");
    fprintf(fp, "# Sixteen rows per read, columns read then 12 of the possible previous base+substitution+next base combinations and the corresponding count\n");
    for(read=0;read<N_READS;read++)
    {
        /* read summary */
        SurvTable *st = read_sts[read];
        if( NULL == st) continue;
        int cntxt_off = 0;
        int len_cntxt = 4;
        int num_cntxt = 256;
        long *count = (long *)smalloc(num_cntxt * sizeof(long));
        for (j=0;j<num_cntxt;j++)
            count[j]=0;
        for (j=0;j<NUM_CNTXT;j++)
        {
            char *cntxt;
            int word;
            cntxt=word2str(j,LEN_CNTXT);
            cntxt+=cntxt_off;
            word=str2word(cntxt,len_cntxt);
            count[word] += st->cntxtH[j];
        }
        for (j=0,p=0;j<num_cntxt;j++)
        {
            char *cntxt;
            cntxt=word2str(j,len_cntxt);
            if (p != cntxt[1]) {
                p = cntxt[1];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "PRCNH\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, count[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }
    
    fprintf(fp, "# Effect of previous base and next base low quality. Use `grep ^PRCNL | cut -f 2-` to extract this part\n");
    fprintf(fp, "# Sixteen rows per read, columns read then 12 of the possible previous base+substitution+next base combinations and the corresponding count\n");
    for(read=0;read<N_READS;read++)
    {
        /* read summary */
        SurvTable *st = read_sts[read];
        if( NULL == st ) continue;
        int cntxt_off = 0;
        int len_cntxt = 4;
        int num_cntxt = 256;
        long *count = (long *)smalloc(num_cntxt * sizeof(long));
        for (j=0;j<num_cntxt;j++)
            count[j]=0;
        for (j=0;j<NUM_CNTXT;j++)
        {
            char *cntxt;
            int word;
            cntxt=word2str(j,LEN_CNTXT);
            cntxt+=cntxt_off;
            word=str2word(cntxt,len_cntxt);
            count[word] += st->cntxtL[j];
        }
        for (j=0,p=0;j<num_cntxt;j++)
        {
            char *cntxt;
            cntxt=word2str(j,len_cntxt);
            if (p != cntxt[1]) {
                p = cntxt[1];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "PRCNL\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, count[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }
    
    for (read=0; read < N_READS; read++) {
        freeSurvEntry(read_sts[read]);
        free(read_sts[read]);
    }
    free(read_sts);
    if (opts->reportName) fclose(fp);
}

static int updateSurvTable(opts_t *opts, SurvTable **sts,
                           int read, int *read_mismatch,
                           char *read_seq, int *read_qual, char *read_ref) {

    int read_length = opts->read_length[read];
    int b;

    /* update survival table */
    for (b = 0; b < read_length; b++) {
        SurvTable *st = sts[read] + b;
        float predictor = -1.0;
        int ibin;

        predictor = read_qual[b];
        
        ibin = predictor + 0.5;
        if( ibin >= st->nbins ) ibin=(st->nbins-1);

        if( read_mismatch[b] & BASE_KNOWN_SNP ) {
            // don't count these
        } else {
            if( read_mismatch[b] & BASE_ALIGN )
                st->num_bases[ibin]++;
            if( read_mismatch[b] & BASE_MISMATCH ){
                char subst[LEN_SUBST+1];
                char cntxt[LEN_CNTXT+1];
                int chr, word;

                st->num_errors[ibin]++;

                chr=0;
                subst[chr++]=read_ref[b];
                subst[chr++]=read_seq[b];
                word=str2word(subst, LEN_SUBST);
                if( word >= 0 ) {
                    st->subst[word][ibin]++;
                    if( predictor >= st->predictor_hilo )
                        st->substH[word]++;
                    else
                        st->substL[word]++;
                }

                chr=0;
                cntxt[chr++]=(b > 0 ? read_ref[b-1] : 'N');
                cntxt[chr++]=read_ref[b];
                cntxt[chr++]=read_seq[b];
                cntxt[chr++]=(b < (read_length-1) ? read_ref[b+1] : 'N');
                word=str2word(cntxt, LEN_CNTXT);
                if( word >= 0 ) {
                    if( predictor >= st->predictor_hilo )
                        st->cntxtH[word]++;
                    else
                        st->cntxtL[word]++;
                }
            }
        }
    }
    
    return 0;
}

/*
 * Takes the bam file as input and updates the survival table
 *
 * Assumption: within a single input file, all reads are the same length and
 * we're using unclipped data.
 *
 * Returns: +'ve integer for success
 *          0 for failure
 */
SurvTable **LoadData(opts_t *opts) {
    SurvTable **sts = NULL;

    sts = smalloc(N_READS * sizeof(SurvTable *));
    for(int read=0;read<N_READS;read++) sts[read] = NULL;

    size_t nreads = 0;

    static const int bam_read_buff_size = 1024;
    char bam_read_seq[bam_read_buff_size];
    int bam_read_qual[bam_read_buff_size];
    char bam_read_ref[bam_read_buff_size];
    int bam_read_mismatch[bam_read_buff_size];

    BAMit_t *bam_in = BAMit_open(opts->in_bam_file, 'r', opts->input_fmt, 0, NULL);
    if (NULL == bam_in) {
        die("ERROR: can't open bam file %s: %s\n", opts->in_bam_file, strerror(errno));
    }

    bam1_t *bam;

    /* loop over reads in the bam file */
    while (1) {
        int bam_read = -1, read_length;

        bam = parse_bam_readinfo(bam_in, NULL, NULL, NULL, NULL, &bam_read, NULL);
        if (!bam) break;    // exit loop at end of BAM file

        if (BAM_FUNMAP & bam->core.flag) continue;
        if (BAM_FQCFAIL & bam->core.flag) continue;
        if (BAM_FSECONDARY & bam->core.flag) continue;
        if (BAM_FSUPPLEMENTARY & bam->core.flag) continue;
        if (BAM_FPAIRED & bam->core.flag) {
            if (BAM_FMUNMAP & bam->core.flag) continue;
            if (0 == (BAM_FPROPER_PAIR & bam->core.flag)) continue;
        }

        read_length = bam->core.l_qseq;
        if (0 == opts->read_length[bam_read]) {
            opts->read_length[bam_read] = read_length;
        }

        if (opts->read_length[bam_read] != read_length) {
            fprintf(stderr,
                    "Error: inconsistent read lengths "
                    "within bam file for read %d.\n"
                    "have length %ld, previously it was %d.\n",
                    bam_read, (long) read_length, opts->read_length[bam_read]);
            exit(EXIT_FAILURE);
        }

        parse_bam_alignments(bam_in, bam, bam_read_seq, bam_read_qual, bam_read_ref,
                                     bam_read_mismatch, bam_read_buff_size, NULL);

        if (NULL == sts[bam_read]) {
            sts[bam_read] = (SurvTable *)smalloc(read_length * sizeof(SurvTable));
            for(int cycle=0;cycle<read_length;cycle++)
                initialiseSurvTable(sts[bam_read]+cycle, bam_read, cycle);
        }
        
        if (0 != updateSurvTable(opts, sts,
                                 bam_read, bam_read_mismatch,
                                 bam_read_seq, bam_read_qual, bam_read_ref)) {
            fprintf(stderr,"ERROR: updating survival table for tile %i.\n", 0);
            exit(EXIT_FAILURE);
        }
        
        nreads++;
    }
    
    completeSurvTable(opts, sts);

    bam_destroy1(bam);
    BAMit_free(bam_in);
    return sts;
}

static void usage(FILE *usagefp)
{
    fprintf(usagefp, "Usage: bambi substitution_analysis [options] bam_file\n");
    fprintf(usagefp, "\n");
    fprintf(usagefp, "Reads the given BAM (or SAM or CRAM) file and produces a substitution analysis table\n");
    fprintf(usagefp, "\n");
    fprintf(usagefp, "Options:\n");
    fprintf(usagefp, " -v --verbose     display progress messages to stderr\n");
    fprintf(usagefp, " -o               output filename for report [default: stdout]\n");
    fprintf(usagefp, "    --input-fmt   BAM input format [sam|bam|cram] [default: bam]\n");
}

/*
 * Parse the command line arguments into a form we can use
 */
opts_t* substitution_analysis_parse_args(int argc, char *argv[])
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char *optstring = "vh?o:";

    static const struct option lopts[] = {
        {"help", 0, 0, 'h'},
        {"verbose", 0, 0, 'v'},
        {"output", 1, 0, 'o'},
        {"output-fmt", 1, 0, 0},
        {"input-fmt", 1, 0, 0},
        {"compression-level", 1, 0, 0},
        {0, 0, 0, 0}
    };

    opts_t *opts = smalloc(sizeof(opts_t));

    for(int n=0; n<N_READS; n++) opts->read_length[n] = 0;
    opts->in_bam_file = NULL;
    opts->reportName = NULL;
    opts->input_fmt = NULL;

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, optstring, lopts, &option_index)) != -1) {
        const char *arg;
        switch (opt) {
            case 'o': opts->reportName = strdup(optarg);      break;
            case 'v': opts->verbose = 1;                    break;
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

    if (optind < argc) opts->in_bam_file = strdup(argv[optind]);

    return opts;
}

static int substitution_analysis(opts_t *opts)
{
    SurvTable **sts = NULL;

    sts = LoadData(opts);
    writeReport(opts, sts);
    freeSurvTable(opts, sts);

    return EXIT_SUCCESS;
}

/*
 * Called from bambi to create a substitution analysis table
 *
 * parse the command line arguments, then call the main process
 *
 * returns 0 on success, 1 if there was a problem
 */
int main_substitution_analysis(int argc, char *argv[])
{
    int ret = 1;
    opts_t *opts = substitution_analysis_parse_args(argc, argv);
    if (opts) ret = substitution_analysis(opts);
    if (opts) free_opts(opts);
    return ret;
}

