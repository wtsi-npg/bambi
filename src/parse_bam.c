/*  parse_bam.c - BAM parsing functions used by spatial_filter

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
#include <assert.h>
#include <string.h>

#include "parse_bam.h"
#include "bambi.h"

#define bam_nt16_rev_table "=ACMGRSVTWYHKDBN"

static char *complement_table = NULL;

/**
 * Reverses the direction of the array of ints
 *
 * @param int is the input array. <b>NB.</b> this is destructively modified.
 *
 * @returns a pointer to the original storage but with the contents
 * modified
 */
int *reverse_int(int *num, int n)
{
    int *t = num, *s = num + n - 1;
    int c;

    while (t < s) {
        c = *s;
        *s = *t;
        *t = c;
        t++;
        s--;
    }

    return num;
}

/**
 * Reverses the direction of the string.
 *
 * @param seq is the input sequence. <b>NB.</b> this is destructively modified.
 *
 * @returns a pointer to the original string storage but with the contents
 * modified
 */
char *reverse_seq(char *seq)
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

/**
 * Return a character representing the complement-base of the supplied
 * parameter.  The mapping is as follows: a->t, c->g, g->c, t|u->a, [->],
 * ]->[, -->-, all others->n. There is a single shot initalisation
 * of a static lookup table to represent this mapping. In addition the
 * case of the supplied parameter is preserved, Uppercase->Uppercase and
 * vice-versa.
 *
 * @param c is the character representing a base. Mapped values include:
 * {a,c,g,t,u,[,],-} all other inputs get a default mapping of 'n'.
 *
 * @returns the character representation of the bilogical compliment of the
 * supplied base. 
 */
char complement_base(char c)
{
    if (!complement_table) {
        int x;
        complement_table = (char *) calloc(256, sizeof(char)) + 127;

        for (x = -127; x < 128; x++) {
                 if (x == 'a') complement_table[x] = 't';
            else if (x == 'c') complement_table[x] = 'g';
            else if (x == 'g') complement_table[x] = 'c';
            else if (x == 't') complement_table[x] = 'a';
            else if (x == 'u') complement_table[x] = 'a';
            else if (x == 'n') complement_table[x] = 'n';
            else if (x == 'A') complement_table[x] = 'T';
            else if (x == 'C') complement_table[x] = 'G';
            else if (x == 'G') complement_table[x] = 'C';
            else if (x == 'T') complement_table[x] = 'A';
            else if (x == 'U') complement_table[x] = 'A';
            else if (x == 'N') complement_table[x] = 'N';
            else complement_table[x] = x;
        }
    }

    return complement_table[(int) c];
}

/**
 * Convert a string containing a string representation of a base sequence
 * into its biological complement. See complement_base for the mapping. This
 * does not reverses the direction of the string.
 *
 * @see complement_base
 *
 * @param seq is the input sequence. <b>NB.</b> this is destructively modified.
 *
 * @returns a pointer to the original string storage but with the contents
 * modified to have complement characters.
 */
char *complement_seq(char *seq)
{
    char *s = seq;
    while (*s) {
        *s = complement_base(*s);
        s++;
    }

    return seq;
}


/* cts -simplification of parse_4_int code for single int parse
 */
const char *parse_next_int(const char *str, int *val, const char *sep)
{
    const static char spaces[256] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    const char* const start = str;
    int minus = 0;
    int ival = 0;
    char c;

    if (NULL == sep) {
        while (*str && spaces[(unsigned) *str]) ++str;
    } else {
        while (*str && NULL != strchr(sep, *str)) ++str;
    }

    c = *str;

    if (!c) {
      /*
        fprintf(stderr, "Error: expected to parse int from string \"%s\"\n", start);
        exit(EXIT_FAILURE);
    */
      return NULL;
    }

    if (c == '-' || c == '+') {
        minus = (c == '-');
        c = *++str;
    }

    while (c >= '0' && c <= '9') {
        ival = ival * 10 + (c-'0');
        c = *++str;
    }

    if (NULL == sep) {
        switch(c) {
        case '\n': case '\r': case '\t': case ' ': case '\0':
            if (minus) ival = -ival;
            *val = ival;
            return str;
        }
    } else {
        if (NULL != strchr(sep, *str)) {
            if (minus) ival = -ival;
            *val = ival;
            return str;
        }
    }
    fprintf(stderr, "Error: expected to parse int from string \"%s\"\n", start);
    exit(EXIT_FAILURE);
}

/*
 * cts - parse bam file line
 *
 * returns 0 on success, 1 on expected failure.
 */
bam1_t *parse_bam_readinfo( BAMit_t *fp,
                    int *bam_lane,
                    int *bam_tile,
                    int *bam_x,
                    int *bam_y,
                    int *bam_read,
                    size_t *bam_offset) 
{

    const char* const sep = ":#/";
    uint8_t *ci_ptr;
    int lane, tile, x, y, read;
	int offset;

    bam1_t *bam = BAMit_next(fp);
    if (!bam) return bam;

    lane = -1;
    tile = -1;
    x = -1;
    y = -1;

    const char* const name = bam_get_qname(bam);
    const char *cp = name;
    /* let's find the beginning of the last 4 subfields of name separated by ':' */
    const char *c_subfield_p[4] = {name, name, name, name};
    int c_subfield_i = 0;
    while (NULL != (cp = strchr(cp,':'))) {
      c_subfield_p[c_subfield_i] = ++cp;
      c_subfield_i++;
      c_subfield_i %= 4;
    }
    cp = c_subfield_p[c_subfield_i];
    cp = parse_next_int(cp,&lane,sep);
    cp = parse_next_int(cp,&tile,sep);
    cp = parse_next_int(cp,&x,sep);
    cp = parse_next_int(cp,&y,sep);

    if (bam_offset) {	/* look for offset, if we want it */
      /* look for ci tag */
      ci_ptr = bam_aux_get(bam, "ci");
      if (NULL == ci_ptr) {
        /* no ci tag get offset from name (cannot be spearated with a ':')*/
        cp = parse_next_int(cp,&offset,sep);
        if (NULL == cp) die("ERROR: No ci tag and no offset in name: \"%s\"\n",name);
      } else {
        offset = bam_aux2i(ci_ptr);
        /* name offset is 0 based but ci is 1 based */
        offset--;
      }
    }

    if (lane < 1) die("ERROR: Invalid lane value in name: \"%s\"\n",name);

    if (tile <= 0) die("ERROR: Invalid tile value in name: \"%s\"\n",name);

    read = 0;
    if (BAM_FPAIRED & bam->core.flag){
        if (BAM_FREAD1 & bam->core.flag) read = 1;
        if (BAM_FREAD2 & bam->core.flag) read = 2;
        if (read == 0){
            die("ERROR: Unable to determine read from flag %d for read: \"%s\"\n",bam->core.flag,name);
        }
    }

    *bam_lane = lane;
    *bam_tile = tile;
    *bam_x = x;
    *bam_y = y;
    *bam_read = read;
    if (bam_offset) *bam_offset = offset;
    return bam;
}


/*
 * cts - parse bam file line
 *
 * returns 0 on success, 1 on expected failure
 */
int parse_bam_alignments(
             BAMit_t *fp,
             bam1_t *bam,
             char *read_seq,
             int *read_qual,
             char *read_ref,
             int *read_mismatch, const int read_buff_size,
	HashTable *snp_hash)
{

    char *name;
    int32_t pos;
    uint32_t *cigar;
    uint8_t *seq, *qual, *m_ptr;
    char *mismatch = NULL;
    const char *sep = "^ACGTKMRYSWBVHDNacgtkmryswbvhdn";
    const char *cp;
    int i, j, skip;
    HashItem *hi;

    if (0 == read_buff_size) die("ERROR: Invalid read_buff_size");

    name = bam_get_qname(bam);
    pos = bam->core.pos;
    cigar = bam_get_cigar(bam);
    seq = bam_get_seq(bam);
    qual = bam_get_qual(bam);
    m_ptr = bam_aux_get(bam, "MD");
    assert(read_buff_size > bam->core.l_qseq);

    if (NULL == m_ptr) {
        die("ERROR: No mismatch for read: \"%s\"\n", name);
    } else {
        mismatch = bam_aux2Z(m_ptr);
        if (NULL == mismatch) {
            die("ERROR: Invalid mismatch %s for read: \"%s\"\n", mismatch, name);
        }
    }

    memset(read_mismatch, 0, bam->core.l_qseq * sizeof(int));

    for (i = 0; i < bam->core.l_qseq; i++) {
        read_seq[i] = bam_nt16_rev_table[bam_seqi(seq, i)];
        read_qual[i] = qual[i];
        if (read_ref) read_ref[i] = 'N';
    }
    read_seq[i] = 0;
    if (read_ref) read_ref[i] = 0;

    j = 0;
    for (i = 0; i < bam->core.n_cigar; i++) {
        int l = cigar[i] >> 4, op = cigar[i] & 0xf, jdel, k;
        switch (op) {
        case BAM_CMATCH:
            // CIGAR: alignment match;
            for (k = 0; k < l; j++, k++) {
                read_mismatch[j] |= BASE_ALIGN;
                if (read_ref) read_ref[j] = read_seq[j];
			}
            break;
        case BAM_CDEL:
            // CIGAR: deletion from the reference
            // deleted bases are missing so tag the cycle immediately before the deletion, since the cigar is in the
            // forward direction the cycle before the deletion is j-1 for forward reads and j for reverse reads
            jdel = j - (BAM_FREVERSE & bam->core.flag ? 0 : 1);
            if (jdel < 0 || jdel >= bam->core.l_qseq) {
                die("ERROR: Deletion at start/end of read: %s\n", name);
            }
            read_mismatch[jdel] |= BASE_DELETION;
            break;
        case BAM_CINS:
            // CIGAR: insertion to the reference 
            for (k = 0; k < l; j++, k++) {
                read_mismatch[j] |= BASE_INSERTION;
				if (read_ref) read_ref[j] = 'I';
			}
            break;
        case BAM_CSOFT_CLIP:
            // CIGAR: clip on the read with clipped sequence present in qseq
            for (k = 0; k < l; j++, k++) {
                read_mismatch[j] |= BASE_SOFT_CLIP;
				if (read_ref) read_ref[j] = 'S';
			}
            break;
        case BAM_CREF_SKIP:
            // CIGAR: skip on the reference (e.g. spliced alignment)
            break;
        default:
            die("ERROR: Unexpected CIGAR operation: %d\n", op);
        }
    }
    if (j != bam->core.l_qseq) {
        die("ERROR: Inconsistent cigar string %d > %d for read: \"%s\"\n", j, bam->core.l_qseq, name);
    }

    /* clipped sequence or insertions are missing from MD */
    for (i = 0, skip = 0; i < bam->core.l_qseq; i++, skip++)
      if (0 == (read_mismatch[i] & (BASE_SOFT_CLIP | BASE_INSERTION)))
            break;

    cp = mismatch;
    while (NULL != (cp = parse_next_int(cp, &j, sep))) {
        /* skip matching bases, exclude insertions which are missing from MD */
        for (; j > 0; i++)
            if (read_mismatch[i] & BASE_INSERTION)
                skip++;
            else
                j--;

        if (0 == strlen(cp))
            /* reached end of MD string */
            break;

        /* skip insertions which are missing from MD */
        for (; i < bam->core.l_qseq; i++, skip++)
            if (0 == (read_mismatch[i] & BASE_INSERTION))
                break;
        if (i == bam->core.l_qseq) {
            die("ERROR: Invalid MD string %s for read: \"%s\"\n",
                mismatch, name);
        }

        switch (*cp) {
        case '^':
            /* deletions missing from read_seq */
            while(NULL != strchr(sep, *(++cp)))
                skip--;
            break;
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            /* mismatch */
            if (0 == (read_mismatch[i] & BASE_ALIGN)) {
                die("ERROR: Inconsistent cigar string expect alignment at mismatch for read: \"%s\"\n", name);
            }
			if (read_ref) read_ref[i] = *cp;
            read_mismatch[i] |= BASE_MISMATCH;

            pos = bam->core.pos + (i - skip);

            if (NULL != snp_hash) {
                char *chrom = fp->h->target_name[bam->core.tid];
                char key[100];
                /* N.B bam->core.pos is 0 based */
                snprintf(key, sizeof(key), "%s:%d", chrom, pos);
                if (NULL !=
                    (hi = HashTableSearch(snp_hash, key, strlen(key)))) {
                    hi->data.i++;
                    read_mismatch[i] |= BASE_KNOWN_SNP;
                }
            }
            i++;
            break;
        default:
            /* treat all other reference bases as a known SNP */
            while(NULL != strchr(sep, *cp)) {
                read_mismatch[i++] |= BASE_KNOWN_SNP;
                cp++;
            }
            break;
        }
    }

    /* clipped sequence or insertions are missing from MD */
    for (; i < bam->core.l_qseq; i++)
        if (0 == (read_mismatch[i] & (BASE_SOFT_CLIP | BASE_INSERTION)))
            break;
    if (i != bam->core.l_qseq) {
        die("ERROR: Inconsistent MD string %d != %d for read: \"%s\"\n",
            i, bam->core.l_qseq, name);
    }

    if (BAM_FREVERSE & bam->core.flag) {
        read_seq = reverse_seq(read_seq);
        read_seq = complement_seq(read_seq);
        read_qual = reverse_int(read_qual, bam->core.l_qseq);
		if (read_ref) {
			read_ref = reverse_seq(read_ref);
			read_ref = complement_seq(read_ref);
		}
        read_mismatch = reverse_int(read_mismatch, bam->core.l_qseq);
    }

    return 0;
}

/*
 * bam_header_add_pg - add PG record to BAM header
 */

static void append_header_text(bam_hdr_t *header, char* text, int len)
{
	int x = header->l_text + 1;
	int y = header->l_text + len + 1; // 1 byte null
	if (text == 0) return;
	if (x < y) header->text = (char*)realloc(header->text, y);
	strncpy(header->text + header->l_text, text, len); // we cannot use strcpy() here.
	header->l_text += len;
	header->text[header->l_text] = 0;
}

void bam_header_add_pg(char *id, char *pn, char *ds, char *cl, bam_hdr_t *bam_header)
{
    char *text, *id2, *pp = NULL, *pg;
    int id2size, id2num, pgsize, pglen;

	if (NULL == bam_header) {
		die("ERROR: No bam header\n");
	}

	if (NULL == bam_header->text) {
		die("ERROR: No text in bam header\n");
	}
	text = strdup(bam_header->text);

        id2size = 10 + strlen(id);
        id2 = malloc(id2size);
        id2num = 0;
        while(id2num < 10) {
                char *hl, *endl, *endt;
                strcpy(text, bam_header->text);
                sprintf(id2, (id2num > 0 ? "%s%d" : "%s"), id, id2num);
                hl = text;
                while (0 < strlen(hl)) {
                        if (NULL == (endl = strchr(hl, '\n'))) {
                                die("ERROR: Corrupt bam header \"%s\"\n", hl);
                        }
                        *endl = 0;
                        if (0 == memcmp(hl, "@PG", 3)) {
                                if (NULL == (pp = strstr(hl, "ID:"))) {
                                        die("ERROR: No ID in PG line \"%s\"\n", hl);
                                }
                                pp += 3;
                                if (NULL != (endt = strchr(pp, '\t'))) {
                                        *endt = 0;
                                }
                                if (0 == strcmp(pp, id2)) {
                                        /* an existing ID matches the new ID */
                                        break;
                                }
                        }
                        hl = endl + 1;
                }
                if( 0 == strlen(hl)) {
                        /* the new ID doesn't match an existing ID */
                        break;
                }
                /* the new ID matches an existing ID, increment the number and try again */
                id2num++;
        }
        if (id2num == 10) {
            die("ERROR: Header already contains PG lines with ID=%s .. ID=%s\n", id, id2);
        }

	pgsize = 128 + strlen(pn) + strlen(pn) + strlen(ds) + strlen(bambi_version()) + strlen(cl);
	if (NULL != pp) {
		pgsize += strlen(pp);
		pg = malloc(pgsize);
		pglen =
		    snprintf(pg, pgsize,
			     "@PG\tID:%s\tPN:%s\tPP:%s\tDS:%s\tVN:%s\tCL:%s\n", id2, pn, pp, ds, bambi_version(), cl);
	} else {
		pg = malloc(pgsize);
		pglen =
		    snprintf(pg, pgsize,
			     "@PG\tID:%s\tPN:%s\tDS:%s\tVN:%s\tCL:%s\n", id2, pn, ds, bambi_version(), cl);
	}
	assert(pglen < pgsize);

	append_header_text(bam_header, pg, pglen);

	free(pg);
        free(id2);
        free(text);
}
