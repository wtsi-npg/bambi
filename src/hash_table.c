/*
 * Copyright (c) 2005-2011, 2017 Genome Research Ltd.
 * Author(s): James Bonfield
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 * 
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *    Institute nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific
 *    prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash_table.h"

/* =========================================================================
 * TCL's hash function. Basically hash*9 + char.
 * =========================================================================
 */

uint32_t HashTcl(uint8_t *data, int len) {
    uint32_t hash = 0;
    int i;

    for (i = 0; i < len; i++) {
	hash += (hash<<3) + data[i];
    }

    return hash;
}

/* =========================================================================
 * Paul Hsieh's hash function
 * http://www.azillionmonkeys.com/qed/hash.html
 * =========================================================================
 */

#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((const uint8_t *)(d))[1] << 8UL)\
                      +((const uint8_t *)(d))[0])
#endif

uint32_t HashHsieh(uint8_t *data, int len) {
    uint32_t hash = 0, tmp;
    int rem;

    if (len <= 0 || data == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
    case 3: hash += get16bits (data);
	hash ^= hash << 16;
	hash ^= data[sizeof (uint16_t)] << 18;
	hash += hash >> 11;
	break;
    case 2: hash += get16bits (data);
	hash ^= hash << 11;
	hash += hash >> 17;
	break;
    case 1: hash += *data;
	hash ^= hash << 10;
	hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 2;
    hash += hash >> 15;
    hash ^= hash << 10;

    return hash;
}

/* =========================================================================
 * Bob Jenkins' hash function
 * http://burtleburtle.net/bob/hash/doobs.html
 *
 * See jenkins_lookup3.c for a new version of this that has good hash
 * characteristics for a full 64-bit hash value.
 * =========================================================================
 */

#define hashsize(n) ((uint32_t)1<<(n))
#define hashmask(n) (hashsize(n)-1)

/*
--------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.
For every delta with one or two bits set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
* If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
* If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
mix() was built out of 36 single-cycle latency instructions in a 
  structure that could supported 2x parallelism, like so:
      a -= b; 
      a -= c; x = (c>>13);
      b -= c; a ^= x;
      b -= a; x = (a<<8);
      c -= a; b ^= x;
      c -= b; x = (b>>13);
      ...
  Unfortunately, superscalar Pentiums and Sparcs can't take advantage 
  of that parallelism.  They've also turned some of those single-cycle
  latency instructions into multi-cycle latency instructions.  Still,
  this is the fastest good hash I could find.  There were about 2^^68
  to choose from.  I only looked at a billion or so.
--------------------------------------------------------------------
*/
#define mix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

/*
--------------------------------------------------------------------
hash() -- hash a variable-length key into a 32-bit value
  k       : the key (the unaligned variable-length array of bytes)
  len     : the length of the key, counting by bytes
  initval : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Every 1-bit and 2-bit delta achieves avalanche.
About 6*len+35 instructions.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (uint8_t **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);

By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.  It's free.

See http://burtleburtle.net/bob/hash/evahash.html
Use for hash table lookup, or anything where one collision in 2^^32 is
acceptable.  Do NOT use for cryptographic purposes.
--------------------------------------------------------------------
*/

uint32_t HashJenkins(uint8_t *k, int length /*, uint32_t initval */)
{
   register uint32_t a,b,c,len;

   /* Set up the internal state */
   len = length;
   a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
   c = 0; /* initval; */        /* the previous hash value */

   /*---------------------------------------- handle most of the key */
   while (len >= 12)
   {
      a += (k[0] +((uint32_t)k[1]<<8) +((uint32_t)k[2]<<16) +((uint32_t)k[3]<<24));
      b += (k[4] +((uint32_t)k[5]<<8) +((uint32_t)k[6]<<16) +((uint32_t)k[7]<<24));
      c += (k[8] +((uint32_t)k[9]<<8) +((uint32_t)k[10]<<16)+((uint32_t)k[11]<<24));
      mix(a,b,c);
      k += 12; len -= 12;
   }

   /*------------------------------------- handle the last 11 bytes */
   c += length;
   switch(len)              /* all the case statements fall through */
   {
   case 11: c+=((uint32_t)k[10]<<24);
   case 10: c+=((uint32_t)k[9]<<16);
   case 9 : c+=((uint32_t)k[8]<<8);
      /* the first byte of c is reserved for the length */
   case 8 : b+=((uint32_t)k[7]<<24);
   case 7 : b+=((uint32_t)k[6]<<16);
   case 6 : b+=((uint32_t)k[5]<<8);
   case 5 : b+=k[4];
   case 4 : a+=((uint32_t)k[3]<<24);
   case 3 : a+=((uint32_t)k[2]<<16);
   case 2 : a+=((uint32_t)k[1]<<8);
   case 1 : a+=k[0];
     /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return c;
}

/*
 * An interface to the above hash functions.
 * Returns:
 *    A 32-bit hash key, suitable for masking down to smaller bit sizes
 */
uint32_t hash(int func, uint8_t *key, int key_len) {
    switch (func) {
    case HASH_FUNC_HSIEH:
	return HashHsieh(key, key_len);

    case HASH_FUNC_TCL:
	return HashTcl(key, key_len);
	
    case HASH_FUNC_JENKINS:
	return HashJenkins(key, key_len);
    }
    
    return 0;
}

/*
 * As per hash() above but returns a 64-bit key. For 32-bit hash functions
 * this is simply a duplication of the 32-bit value.
 */
uint64_t hash64(int func, uint8_t *key, int key_len) {
    uint32_t pc = 0, pb = 0;

    switch (func) {
    case HASH_FUNC_HSIEH:
	pb = pc = HashHsieh(key, key_len);
	break;

    case HASH_FUNC_TCL:
	pb = pc = HashTcl(key, key_len);
	break;
	
    case HASH_FUNC_JENKINS:
	pb = pc = HashJenkins(key, key_len);
	break;
    }
    
    return pc + (((uint64_t)pb)<<32);
}

/* =========================================================================
 * Hash Table handling code
 * =========================================================================
 */

/* Multiplicative factors indicating when to grow or shrink the hash table */
#define HASH_TABLE_RESIZE 3

/*
 * Creates a HashItem for use with HashTable h.
 *
 * Returns:
 *    A pointer to new HashItem on success
 *    NULL on failure.
 */
static HashItem *HashItemCreate(HashTable *h) {
    HashItem *hi;

    hi = malloc(sizeof(*hi));
    if (NULL == hi) return NULL;

    hi->data.p    = NULL;
    hi->data.i    = 0;
    hi->next      = NULL;
    hi->key       = NULL;
    hi->key_len   = 0;

    h->nused++;
    
    return hi;
}

/*
 * Deallocates a HashItem created via HashItemCreate.
 *
 * This function will not remove the item from the HashTable so be sure to
 * call HashTableDel() first if appropriate.
 */
static void HashItemDestroy(HashTable *h, HashItem *hi, int deallocate_data) {
    if (!hi) return;

    if (!(h->options & HASH_NONVOLATILE_KEYS) || (h->options & HASH_OWN_KEYS))
	if (hi->key)
	    free(hi->key);

    if (deallocate_data && hi->data.p) {
        free(hi->data.p);
    }

	free(hi);

    h->nused--;
}

/*
 * Creates a new HashTable object. Size will be rounded up to the next
 * power of 2. It is a starting point and hash tables may be grown or shrunk
 * as needed (if HASH_DYNAMIC_SIZE is used).
 *
 * If HASH_POOL_ITEMS is used, HashItems will be allocated in blocks to reduce
 * malloc overhead in the case where a large number of items is required.
 * HashItems allocated this way will be put on a free list when destroyed; the
 * memory will only be reclaimed when the entire hash table is destroyed.
 *
 * Options are as defined in the header file (see HASH_* macros).
 *
 * Returns:
 *    A pointer to a HashTable on success
 *    NULL on failure
 */
HashTable *HashTableCreate(int size, int options) {
    HashTable *h;
    int i, bits;
    uint32_t mask;

    if (!(h = (HashTable *)malloc(sizeof(*h))))
	return NULL;

    if (size < 4)
	size = 4; /* an inconsequential minimum size */

    /* Round the requested size to the next power of 2 */
    bits = 0;
    size--;
    while (size) {
	size /= 2;
	bits++;
    }
    size = 1<<bits;
    mask = size-1;

    h->nbuckets = size;
    h->mask = mask;
    h->options = options;
    h->nused = 0;
    h->bucket = (HashItem **)malloc(sizeof(*h->bucket) * size);
    if (NULL == h->bucket) {
        HashTableDestroy(h, 0);
        return NULL;
    }

    for (i = 0; i < size; i++) {
	h->bucket[i] = NULL;
    }

    return h;
}

/*
 * Deallocates a HashTable object (created by HashTableCreate).
 *
 * The deallocate_data parameter is a boolean to indicate whether the
 * data attached to the hash table should also be free()d. DO NOT USE
 * this if the HashData attached was not a pointer allocated using
 * malloc().
 */
void HashTableDestroy(HashTable *h, int deallocate_data) {
    int i;

    if (!h)
	return;

    if (h->bucket) {
        for (i = 0; i < h->nbuckets; i++) {
	    HashItem *hi = h->bucket[i], *next = NULL;
	    for (hi = h->bucket[i]; hi; hi = next) {
	        next = hi->next;
		HashItemDestroy(h, hi, deallocate_data);
	    }
	}

	free(h->bucket);
    }

    free(h);
}

/*
 * Resizes a HashTable to have 'newsize' buckets.
 * This is called automatically when adding or removing items so that the
 * hash table keeps at a sensible scale.
 *
 * FIXME: Halving the size of the hash table is simply a matter of coaelescing
 * every other bucket. Instead we currently rehash (which is slower).
 * Doubling the size of the hash table currently requires rehashing, but this
 * too could be optimised by storing the full 32-bit hash of the key along
 * with the key itself. This then means that it's just a matter of seeing what
 * the next significant bit is. It's a memory vs speed tradeoff though and
 * re-hashing is pretty quick.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int HashTableResize(HashTable *h, int newsize) {
    HashTable *h2;
    int i;

    /* fprintf(stderr, "Resizing to %d\n", newsize); */

    /* Create a new hash table and rehash everything into it */
    h2 = HashTableCreate(newsize, h->options);

    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi, *next;
	for (hi = h->bucket[i]; hi; hi = next) {
	    uint64_t hv = h2->options & HASH_INT_KEYS
		? hash64(h2->options & HASH_FUNC_MASK,
			 (uint8_t *)&hi->key, hi->key_len) & h2->mask
		: hash64(h2->options & HASH_FUNC_MASK,
			 (uint8_t *)hi->key, hi->key_len) & h2->mask;
	    next = hi->next;
	    hi->next = h2->bucket[hv];
	    h2->bucket[hv] = hi;
	}
    }

    /* Swap the links over & free */
    free(h->bucket);
    h->bucket   = h2->bucket;
    h->nbuckets = h2->nbuckets;
    h->mask     = h2->mask;

    free(h2);

    return 0;
}

/*
 * Adds a HashData item to HashTable h with a specific key. Key can be binary
 * data, but if key_len is passed as zero then strlen() will be used to
 * determine the key length.
 *
 * The "new" pointer may be passed as NULL. When not NULL it is filled out
 * as a boolean to indicate whether the key is already in this hash table.
 *
 * The HASH_ALLOW_DUP_KEYS option (specified when using HashTableCreate)
 * will allow duplicate keys to be stored, and hence *new is also zero.
 * By default duplicate keys are disallowed.
 *
 * Keys are considered to be volatile memory (ie temporary storage) and so the
 * hash table takes separate copies of them. To avoid this use the
 * HASH_NONVOLATILE_KEYS option.
 *
 * If the HASH_OWN_KEYS option was specified when creating the table then
 * keys will be considered to be owned by the hash table. In this case
 * the key will be freed when the table is destroyed regardless of
 * whether the HASH_NONVOLATILE_KEYS option was used to allocate its
 * own private copy.
 *
 * Returns:
 *    The HashItem created (or matching if a duplicate) on success
 *    NULL on failure
 */
HashItem *HashTableAdd(HashTable *h, char *key, int key_len, HashData data,
		       int *new) {
    uint64_t hv;
    HashItem *hi;

    if (!key_len)
	key_len = strlen(key);

    hv = h->options & HASH_INT_KEYS
	? hash64(h->options & HASH_FUNC_MASK, (uint8_t *)&key, key_len) & h->mask
	: hash64(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len) & h->mask;

    /* Already exists? */
    if (!(h->options & HASH_ALLOW_DUP_KEYS)) {
	for (hi = h->bucket[hv]; hi; hi = hi->next) {
	    if (h->options & HASH_INT_KEYS) {
		if ((int)(size_t)hi->key == (int)(size_t)key) {
		    if (new) *new = 0;
		    return hi;
		}
	    } else {
		if (key_len == hi->key_len && key[0] == hi->key[0] &&
		    memcmp(key, hi->key, key_len) == 0) {
		    if (new) *new = 0;
		    return hi;
		}
	    }
	}
    }

    /* No, so create a new one and link it in */
    if (NULL == (hi = HashItemCreate(h)))
	return NULL;

    if (h->options & HASH_NONVOLATILE_KEYS)
	hi->key = key;
    else {
	hi->key = (char *)malloc(key_len+1);
	memcpy(hi->key, key, key_len);
	hi->key[key_len] = 0; /* null terminate incase others print keys */
    }
    hi->key_len = key_len;
    hi->data = data;
    hi->next = h->bucket[hv];
    h->bucket[hv] = hi;

    if ((h->options & HASH_DYNAMIC_SIZE) &&
	h->nused > HASH_TABLE_RESIZE * h->nbuckets)
	HashTableResize(h, h->nbuckets*4);

    if (new) *new = 1;

    return hi;
}


/*
 * Removes a specified HashItem from the HashTable. (To perform this it needs
 * to rehash based on the hash key as hash_item only has a next pointer and
 * not a previous pointer.)
 * 
 * The HashItem itself is also destroyed (by an internal call to
 * HashItemDestroy). The deallocate_data parameter controls whether the data
 * associated with the HashItem should also be free()d.
 *
 * See also the HashTableRemove() function to remove by key instead of
 * HashItem.
 *
 * Returns 0 on success
 *        -1 on failure (eg HashItem not in the HashTable);
 */
int HashTableDel(HashTable *h, HashItem *hi, int deallocate_data) {
    uint64_t hv;
    HashItem *next, *last;

    hv = h->options & HASH_INT_KEYS
	? hash64(h->options & HASH_FUNC_MASK,
		 (uint8_t *)&hi->key, hi->key_len) & h->mask
	: hash64(h->options & HASH_FUNC_MASK,
		 (uint8_t *)hi->key, hi->key_len) & h->mask;

    for (last = NULL, next = h->bucket[hv]; next;
	 last = next, next = next->next) {
	if (next == hi) {
	    /* Link last to next->next */
	    if (last)
		last->next = next->next;
	    else
		h->bucket[hv] = next->next;

	    HashItemDestroy(h, hi, deallocate_data);

	    return 0;
	}
    }

    return -1;
}


/*
 * Searches the HashTable for the data registered with 'key' and removes
 * these items from the HashTable. In essence this is a combination of
 * HashTableSearch and HashTableDel functions.
 *
 * If HASH_ALLOW_DUP_KEYS is used this will remove all items matching 'key',
 * otherwise just a single item will be removed.
 *
 * If 'deallocate_data' is true the data associated with the HashItem will
 * be free()d.
 *
 * Returns
 *    0 on success (at least one item found)
 *   -1 on failure (no items found).
 */
int HashTableRemove(HashTable *h, char *key, int key_len,
		    int deallocate_data) {
    uint64_t hv;
    HashItem *last, *next, *hi;
    int retval = -1;

    if (!key_len)
	key_len = strlen(key);

    hv = h->options & HASH_INT_KEYS
	? hash64(h->options & HASH_FUNC_MASK, (uint8_t *)&key, key_len) & h->mask
	: hash64(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len) & h->mask;

    last = NULL;
    next = h->bucket[hv];

    while (next) {
	hi = next;
	if (((h->options & HASH_INT_KEYS)
	     ? ((int)(size_t)key == (int)(size_t)hi->key)
	     : (key_len == hi->key_len && memcmp(key, hi->key, key_len) == 0))) {
	    /* An item to remove, adjust links and destroy */
	    if (last)
		last->next = hi->next;
	    else
		h->bucket[hv] = hi->next;

	    next = hi->next;
	    HashItemDestroy(h, hi, deallocate_data);

	    retval = 0;
	    if (!(h->options & HASH_ALLOW_DUP_KEYS))
		break;

	} else {
	    /* We only update last when it's something we haven't destroyed */
	    last = hi;
	    next = hi->next;
	}
    }

    return retval;
}

/*
 * Searches the HashTable for the data registered with 'key'.
 * If HASH_ALLOW_DUP_KEYS is used this will just be the first one found.
 * You will then need to use HashTableNext to iterate through the matches.
 *
 * Returns
 *    HashItem if found
 *    NULL if not found
 */
HashItem *HashTableSearch(HashTable *h, char *key, int key_len) {
    uint64_t hv;
    HashItem *hi;

    if (!key_len)
	key_len = strlen(key);

    if (h->options & HASH_INT_KEYS) {
	hv = hash64(h->options & HASH_FUNC_MASK, (uint8_t *)&key, key_len)& h->mask;

	for (hi = h->bucket[hv]; hi; hi = hi->next) {
	    if ((int)(size_t)key == (int)(size_t)hi->key)
		return hi;
	}
    } else {
	hv = hash64(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len) & h->mask;

	for (hi = h->bucket[hv]; hi; hi = hi->next) {
	    if (key_len == hi->key_len &&
		memcmp(key, hi->key, key_len) == 0)
		return hi;
	}
    }

    return NULL;
}

/*
 * Find the next HashItem (starting from 'hi') to also match this key.
 * This is only valid when the HASH_ALLOW_DUP_KEYS is in use and
 * we're not using HASH_INT_KEYS.
 *
 * Returns
 *    HashItem if found
 *    NULL if not found
 */
HashItem *HashTableNext(HashItem *hi, char *key, int key_len) {
    if (!hi)
	return NULL;

    for (hi = hi->next; hi; hi = hi->next) {
	if (key_len == hi->key_len &&
	    memcmp(key, hi->key, key_len) == 0)
	    return hi;
    }

    return NULL;
}

HashItem *HashTableNextInt(HashItem *hi, char *key, int key_len) {
    if (!hi)
	return NULL;

    for (hi = hi->next; hi; hi = hi->next) {
	if (key_len == hi->key_len &&
	    memcmp(&key, &hi->key, key_len) == 0)
	    return hi;
    }

    return NULL;
}

/*
 * Dumps a textual represenation of the hash table to stdout.
 */
void HashTableDump(HashTable *h, FILE *fp, char *prefix) {
    int i;
    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    if (h->options & HASH_INT_KEYS) {
		fprintf(fp, "%s%d => %"PRId64" (0x%"PRIx64")\n",
			prefix ? prefix : "",
			(int)(size_t)hi->key,
			hi->data.i, hi->data.i);
	    } else {
		fprintf(fp, "%s%.*s => %"PRId64" (0x%"PRIx64")\n",
			prefix ? prefix : "",
			hi->key_len, hi->key,
			hi->data.i, hi->data.i);
	    }
	}
    }
}

/*
 * Produces some simple statistics on the hash table population.
 */
void HashTableStats(HashTable *h, FILE *fp) {
    int i;
    double avg = (double)h->nused / h->nbuckets;
    double var = 0;
    int maxlen = 0;
    int filled = 0;
    int clen[51];

    for (i = 0; i <= 50; i++)
	clen[i] = 0;

    for (i = 0; i < h->nbuckets; i++) {
	int len = 0;
	HashItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    len++;
	}
	if (len > 0) {
	    filled++;
	    if (len > maxlen)
		maxlen = len;
	}
	clen[len <= 50 ? len : 50]++;
	var += (len-avg) * (len-avg);
    }
    var /= h->nbuckets;
    /* sd = sqrt(var); */

    fprintf(fp, "Nbuckets  = %d\n", h->nbuckets);
    fprintf(fp, "Nused     = %d\n", h->nused);
    fprintf(fp, "Avg chain = %f\n", avg);
    fprintf(fp, "Chain var.= %f\n", var);
    fprintf(fp, "%%age full = %f\n", (100.0*filled)/h->nbuckets);
    fprintf(fp, "max len   = %d\n", maxlen);
    for (i = 0; i <= maxlen; i++) {
	fprintf(fp, "Chain %2d   = %d\n", i, clen[i]);
    }
}

/*
 * Iterates through members of a hash table returning items sequentially.
 *
 * Returns the next HashItem on success
 *         NULL on failure.
 */
HashItem *HashTableIterNext(HashTable *h, HashIter *iter) {
    do {
	if (iter->hi == NULL) {
	    if (++iter->bnum >= h->nbuckets)
		break;
	    iter->hi = h->bucket[iter->bnum];
	} else {
	    iter->hi = iter->hi->next;
	}
    } while (!iter->hi);
    
    return iter->hi;
}

void HashTableIterReset(HashIter *iter) {
    if (iter) {
	iter->bnum = -1;
	iter->hi = NULL;
    }
}

HashIter *HashTableIterCreate(void) {
    HashIter *iter = (HashIter *)malloc(sizeof(*iter));

    HashTableIterReset(iter);
    return iter;
}

void HashTableIterDestroy(HashIter *iter) {
    if (iter)
	free(iter);
}
