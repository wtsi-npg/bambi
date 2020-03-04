/*  decode.h -- barcode tag decoding

    Copyright (C) 2018 Genome Research Ltd.

    Author: Rob Davies <rmd@sanger.ac.uk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DECODE_H
#define DECODE_H

#include "array.h"
#include "hash_table.h"
#include <stdbool.h>

// Options
typedef struct decode_opts_t decode_opts_t;

decode_opts_t *decode_init_opts(int argc, char **argv);
void decode_free_opts(decode_opts_t *opts);
void set_decode_opt_barcode_name(decode_opts_t *opts, const char *name);
void set_decode_opt_metrics_name(decode_opts_t *opts, const char *name);
void set_decode_opt_barcode_tag_name(decode_opts_t *opts, const char *name);
void set_decode_opt_max_low_quality_to_convert(decode_opts_t *opts, int val);
void set_decode_opt_convert_low_quality(decode_opts_t *opts, bool flag);
void set_decode_opt_max_no_calls(decode_opts_t *opts, int val);
void set_decode_opt_max_mismatches(decode_opts_t *opts, int val);
void set_decode_opt_min_mismatch_delta(decode_opts_t *opts, int val);
void set_decode_opt_change_read_name(decode_opts_t *opts, bool flag);
void set_decode_opt_ignore_pf(decode_opts_t *opts, bool flag);

// Read the barcode definitions file
va_t *loadBarcodeFile(decode_opts_t *opts);

// Create a hash table to speed up barcode searches
HashTable *make_barcode_hash(va_t *barcodeArray);

// Returns the length of the longest name
size_t find_longest_barcode_name(va_t *barcodeArray);

// Get barcode metadata.  Returns -1 if idx is off the end of barcodeArray, else 0
int get_barcode_metadata(va_t *barcodeArray, int idx,
                         const char **name_out, const char **lib_out,
                         const char **sample_out, const char **desc_out,
                         const char **seq_out);

// Make a copy of a barcode array
va_t *copy_barcode_array(va_t *barcode_array);

// Clean up a barcode array copy
void delete_barcode_array_copy(va_t *barcode_array);

// Barcode look-up.  Also updates metrics.
char *findBarcodeName(char *barcode, va_t *barcodeArray, HashTable *barcodeHash, HashTable *tagHopHash, decode_opts_t *opts, bool isPf, bool isUpdateMetrics);

// Accumulate metrics in per-thread copies
void accumulate_job_metrics(va_t *job_barcodes, HashTable *job_tag_hops, va_t *barcodeArray, HashTable *tagHopHash);

// Write out metrics
int writeMetrics(va_t *barcodeArray, HashTable *tagHopHash, decode_opts_t *opts);

// free the tag hop hash table
void free_tagHopHash(HashTable *tagHopHash);

#endif
