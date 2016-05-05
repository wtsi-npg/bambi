/*  i2b.c -- index i2br subcommand.

    Copyright (C) 2016 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include "bambi.h"
#include <assert.h>
#include <ctype.h>
#include <htslib/sam.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <unistd.h>
#include <regex.h>
#include <libgen.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>

#include <cram/sam_header.h>

#define DEFAULT_BARCODE_TAG "BC"
#define DEFAULT_QUALITY_TAG "QT"

/*
 * index array type
 */
typedef struct {
    int end;
    int max;
    int *entries;
} ia_t;

/*
 * structure to hold options
 */
typedef struct {
    int verbose;
    char *argv_list;
    char *run_folder;
    char *intensity_dir;
    char *basecalls_dir;
    int lane;
    char *output_file;
    char *output_fmt;
    char compression_level;
    bool generate_secondary_basecalls;
    bool no_filter;
    char *read_group_id;
    char *sample_alias;
    char *library_name;
    char *study_name;
    char *platform_unit;
    char *run_start_date;
    char *sequencing_centre;
    char *platform;
    int first_tile;
    int tile_limit;
    char *barcode_tag;
    char *quality_tag;
    char *barcode_tag2;
    char *quality_tag2;
    int bc_read;
    int sec_bc_read;
    ia_t *first_cycle;
    ia_t *final_cycle;
    ia_t *first_index_cycle;
    ia_t *final_index_cycle;
    bool add_cluster_index_tag;
    xmlDocPtr intensityConfig;
    xmlDocPtr basecallsConfig;
    xmlDocPtr parametersConfig;
    xmlDocPtr runinfoConfig;
} opts_t;


//
// TODO
// These array handling functions really need to go into a separate module, with seperate tests
//

/*
 * integer array functions
 */

int ia_compare(const void *ia1, const void *ia2)
{
    return *(int *)ia1 - *(int *)ia2;
}

void ia_sort(ia_t *ia)
{
    qsort(ia->entries, ia->end, sizeof(int), ia_compare);
}

void ia_push(ia_t *ia, int i)
{
    if (ia->end == ia->max) {
        // expand the array
        ia->max *= 2;
        ia->entries = realloc(ia->entries, ia->max * sizeof(int));
    }
    ia->entries[ia->end] = i;
    ia->end++;
}

void ia_free(ia_t *ia)
{
    free(ia->entries);
    free(ia);
}

ia_t *ia_init(int max)
{
    ia_t *ia = calloc(1, sizeof(ia_t));
    ia->end = 0;
    ia->max = max;
    ia->entries = calloc(ia->max, sizeof(int));
    return ia;
}

/*
 * Cycle range array
 */

typedef struct {
    char *readname;
    int first, last;
} cycleRangeEntry_t;

void freeCycleRange(void *ent)
{
    cycleRangeEntry_t *cr = (cycleRangeEntry_t *)ent;
    free(cr->readname);
    free(cr);
}

/*
 * generic arrays
 */
typedef struct {
    int end;
    int max;
    void (*free_entry)(void *);
    void **entries;
} va_t;

va_t *va_init(int max, void(*free_entry)(void*))
{
    va_t *va = calloc(1,sizeof(va_t));
    va->end = 0;
    va->max = max;
    va->free_entry = free_entry;
    va->entries = calloc(va->max, sizeof(void *));
    return va;
}

void va_push(va_t *va, void *ent)
{
    if (va->end == va->max) {
        // expand the array
        va->max *= 2;
        va->entries = realloc(va->entries, va->max * sizeof(void *));
    }
    va->entries[va->end] = ent;
    va->end++;
}

bool va_isEmpty(va_t *va)
{
    return va->end == 0;
}

void va_free(va_t *va)
{
    int n;
    for (n=0; n < va->end; n++) {
        va->free_entry(va->entries[n]);
    }
    free(va->entries);
    free(va);
}




static void free_opts(opts_t* opts)
{
    if (!opts) return;
    free(opts->run_folder);
    free(opts->intensity_dir);
    free(opts->basecalls_dir);
    free(opts->argv_list);
    free(opts->output_file);
    free(opts->output_fmt);
    free(opts->read_group_id);
    free(opts->sample_alias);
    free(opts->library_name);
    free(opts->study_name);
    free(opts->platform_unit);
    free(opts->run_start_date);
    free(opts->sequencing_centre);
    free(opts->platform);
    free(opts->barcode_tag);
    free(opts->quality_tag);
    free(opts->barcode_tag2);
    free(opts->quality_tag2);
    ia_free(opts->first_cycle);
    ia_free(opts->final_cycle);
    ia_free(opts->first_index_cycle);
    ia_free(opts->final_index_cycle);
    xmlFreeDoc(opts->intensityConfig);
    xmlFreeDoc(opts->basecallsConfig);
    xmlFreeDoc(opts->parametersConfig);
    xmlFreeDoc(opts->runinfoConfig);
    free(opts);
}

xmlXPathObjectPtr getnodeset(xmlDocPtr doc, char *xpath)
{
    xmlXPathContextPtr context;
    xmlXPathObjectPtr result;

    context = xmlXPathNewContext(doc);
    if (context == NULL) {
        printf("Error in xmlXPathNewContext\n");
        return NULL;
    }
    result = xmlXPathEvalExpression((xmlChar *)xpath, context);
    xmlXPathFreeContext(context);
    if (result == NULL) {
        printf("Error in xmlXPathEvalExpression\n");
        return NULL;
    }
    if(xmlXPathNodeSetIsEmpty(result->nodesetval)){
        xmlXPathFreeObject(result);
                printf("No result\n");
        return NULL;
    }
    return result;
}

static char *getXMLAttr(xmlDocPtr doc, char *node, char *attr)
{
    char *v = NULL;
    xmlNodeSetPtr nodeset;
    xmlXPathObjectPtr result;

    if (!doc) return "";
    result = getnodeset(doc, node);
    if (result) {
        nodeset = result->nodesetval;
        v = (char *)xmlGetProp(nodeset->nodeTab[0], (xmlChar *)attr);
        xmlXPathFreeObject (result);
    }
    return v;
}

xmlDocPtr loadXML(char *dir, char *fname, int verbose)
{
    xmlDocPtr doc;
    char *tmp = calloc(1, strlen(dir) + strlen(fname) + 2);
    sprintf(tmp, "%s/%s", dir, fname);
    doc = xmlReadFile(tmp, NULL, XML_PARSE_NOWARNING);
    if (!doc) {
        if (verbose) fprintf(stderr, "WARNING: Failed to parse %s/%s\n", dir, fname);
    }
    free(tmp);
    return doc;
}

void dumpOpts(opts_t *o)
{
    printf("OptionsL\n");
    printf("verbose:         %d\n", o->verbose);
    printf("argv_list:       %s\n", o->argv_list);
    printf("run_folder:      %s\n", o->run_folder);
    printf("intensity_dir:   %s\n", o->intensity_dir);
    printf("basecalls_dir:   %s\n", o->basecalls_dir);
    printf("Lane:            %d\n", o->lane);
    printf("output_file:     %s\n", o->output_file);
    printf("\n");

/*
    char *output_fmt;
    char compression_level;
    bool generate_secondary_basecalls;
    bool no_filter;
    char *read_group_id;
    char *sample_alias;
    char *library_name;
    char *study_name;
    char *platform_unit;
    char *run_start_date;
    char *sequencing_centre;
    char *platform;
    int first_tile;
    int tile_limit;
    char *barcode_tag;
    char *quality_tag;
    char *barcode_tag2;
    char *quality_tag2;
    int bc_read;
    int sec_bc_read;
    ia_t *first_cycle;
    ia_t *final_cycle;
    ia_t *first_index_cycle;
    ia_t *final_index_cycle;
    bool add_cluster_index_tag;
*/
}

/*
 * display usage information
 */
static void usage(FILE *write_to)
{
    fprintf(write_to,
"Usage: bambi i2b [options]\n"
"\n"
"Options:\n"
"  -r   --run-folder                    Illumina runfolder directory including runParameters xml file under it.\n"
"                                       [default: two levels up from Intensities directory]\n"
"  -i   --intensity-dir                 Illumina intensities directory including config xml file, and clocs,\n"
"                                       locs or pos files under lane directory. Required\n"
"  -b   --basecalls-dir                 Illumina basecalls directory including config xml file, and filter files,\n"
"                                       bcl, maybe scl files under lane cycle directory\n"
"                                       [default: BaseCalls directory under intensities]\n"
"  -l   --lane                          Lane number. Required\n"
"  -o   --output-file                   Output file name. May be '-' for stdout. Required\n"
"       --generate-secondary-basecalls  Including second base call or not [default: false]\n"
"       --no-filter                     Do not filter cluster [default: false]\n"
"       --read-group-id                 ID used to link RG header record with RG tag in SAM record. [default: '1']\n"
"       --library-name                  The name of the sequenced library. [default: 'unknown']\n"
"       --sample-alias                  The name of the sequenced sample. [default: same as library name]\n"
"       --study-name                    The name of the study. [default: none]\n"
"       --platform-unit                 The platform unit. [default: runfolder name plus lane number]\n"
"       --run-start-date                The start date of the run [default: read from config file]\n"
"       --sequencing-centre             Sequencing Centre. [default: 'SC']\n"
"       --platform                      Sequencing technology used. [default: 'ILLUMINA']\n"
"       --first-tile                    First tile to be processed. This is normally only used for testing and\n"
"                                       debugging. [default: null]\n"
"       --tile-limit                    Number of tiles to process. Normally only used for testing and\n"
"                                       debugging. [default: all tiles]\n"
"       --barcode-tag                   Tag name for barcode sequence. [default: " DEFAULT_BARCODE_TAG "]\n"
"       --quality-tag                   Tag name for barcode quality. [default: " DEFAULT_QUALITY_TAG "]\n"
"       --sec-barcode-tag               Tag name for second barcode sequence. [default: null]\n"
"       --sec-quality-tag               Tag name for second barcode quality. [default: null]\n"
"       --bc-read                       Which read (1 or 2) should the barcode sequence and quality be added to?\n"
"                                       [default: 1]\n"
"       --sec-bc-read                   Which read (1 or 2) should the second barcode sequence and quality be added to?\n"
"                                       [default: bc-read]\n"
"       --first-cycle                   First cycle for each standard (non-index) read. Can be specified 0 or more times.\n"
"       --final-cycle                   Last cycle for each standard (non-index) read. Can be specified 0 or more times.\n"
"       --first-index-cycle             First cycle for each index read. Can be specified 0 or more times.\n"
"       --final-index-cycle             Last cycle for each index read. Can be specified 0 or more times.\n"
"       --add-cluster-index-tag         Add cluster index tag [default: false]\n"
"  -v   --verbose                       verbose output\n"
"       --output-fmt                    [sam/bam/cram] [default: bam]\n"
"       --compression-level             [0..9]\n"
);
}

/*
 * Takes the command line options and turns them into something we can understand
 */
opts_t* i2b_parse_args(int argc, char *argv[])
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char* optstring = "vr:i:b:l:o:";

    static const struct option lopts[] = {
        { "verbose",                    0, 0, 'v' },
        { "run-folder",                 1, 0, 'r' },
        { "intensity-dir",              1, 0, 'i' },
        { "basecalls-dir",              1, 0, 'b' },
        { "lane",                       1, 0, 'l' },
        { "output-file",                1, 0, 'o' },
        { "generate-secondary-basecalls", 0, 0, 0 },
        { "no-filter",                  0, 0, 0 },
        { "read-group-id",              1, 0, 0 },
        { "output-fmt",                 1, 0, 0 },
        { "compression-level",          1, 0, 0 },
        { "library-name",               1, 0, 0 },
        { "sample-alias",               1, 0, 0 },
        { "study-name",                 1, 0, 0 },
        { "platform-unit",              1, 0, 0 },
        { "run-start-date",             1, 0, 0 },
        { "sequencing-centre",          1, 0, 0 },
        { "sequencing-center",          1, 0, 0 },
        { "platform",                   1, 0, 0 },
        { "first-tile",                 1, 0, 0 },
        { "tile-limit",                 1, 0, 0 },
        { "barcode-tag",                1, 0, 0 },
        { "quality-tag",                1, 0, 0 },
        { "sec-barcode-tag",            1, 0, 0 },
        { "sec-quality-tag",            1, 0, 0 },
        { "bc-read",                    1, 0, 0 },
        { "sec-bc-read",                1, 0, 0 },
        { "first-cycle",                1, 0, 0 },
        { "final-cycle",                1, 0, 0 },
        { "first-index-cycle",          1, 0, 0 },
        { "final-index-cycle",          1, 0, 0 },
        { "add-cluster-index-tag",      0, 0, 0 },
        { NULL, 0, NULL, 0 }
    };

    opts_t* opts = calloc(sizeof(opts_t), 1);
    if (!opts) { perror("cannot allocate option parsing memory"); return NULL; }

    opts->argv_list = stringify_argv(argc+1, argv-1);
    if (opts->argv_list[strlen(opts->argv_list)-1] == ' ') opts->argv_list[strlen(opts->argv_list)-1] = 0;

    opts->first_cycle = ia_init(5);
    opts->final_cycle = ia_init(5);
    opts->first_index_cycle = ia_init(5);
    opts->final_index_cycle = ia_init(5);

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, optstring, lopts, &option_index)) != -1) {
        const char *arg;
        switch (opt) {
        case 'r':   opts->run_folder = strdup(optarg);
                    break;
        case 'i':   opts->intensity_dir = strdup(optarg);
                    break;
        case 'b':   opts->basecalls_dir = strdup(optarg);
                    break;
        case 'o':   opts->output_file = strdup(optarg);
                    break;
        case 'l':   opts->lane = atoi(optarg);
                    break;
        case 'v':   opts->verbose++;
                    break;
        case 0:     arg = lopts[option_index].name;
                         if (strcmp(arg, "output-fmt") == 0)                   opts->output_fmt = strdup(optarg);
                    else if (strcmp(arg, "compression-level") == 0)            opts->compression_level = *optarg;
                    else if (strcmp(arg, "generate-secondary-basecalls") == 0) opts->generate_secondary_basecalls = true;
                    else if (strcmp(arg, "no-filter") == 0)                    opts->no_filter = true;
                    else if (strcmp(arg, "read-group-id") == 0)                opts->read_group_id = strdup(optarg);
                    else if (strcmp(arg, "library-name") == 0)                 opts->library_name = strdup(optarg);
                    else if (strcmp(arg, "sample-alias") == 0)                 opts->sample_alias = strdup(optarg);
                    else if (strcmp(arg, "study-name") == 0)                   opts->study_name = strdup(optarg);
                    else if (strcmp(arg, "platform-unit") == 0)                opts->platform_unit = strdup(optarg);
                    else if (strcmp(arg, "run-start-date") == 0)               opts->run_start_date = strdup(optarg);
                    else if (strcmp(arg, "sequencing-centre") == 0)            opts->sequencing_centre = strdup(optarg);
                    else if (strcmp(arg, "sequencing-center") == 0)            opts->sequencing_centre = strdup(optarg);
                    else if (strcmp(arg, "platform") == 0)                     opts->platform = strdup(optarg);
                    else if (strcmp(arg, "first-tile") == 0)                   opts->first_tile = atoi(optarg);
                    else if (strcmp(arg, "tile-limit") == 0)                   opts->tile_limit = atoi(optarg);
                    else if (strcmp(arg, "barcode-tag") == 0)                  opts->barcode_tag = strdup(optarg);
                    else if (strcmp(arg, "quality-tag") == 0)                  opts->quality_tag = strdup(optarg);
                    else if (strcmp(arg, "sec-barcode-tag") == 0)              opts->barcode_tag2 = strdup(optarg);
                    else if (strcmp(arg, "sec-quality-tag") == 0)              opts->quality_tag2 = strdup(optarg);
                    else if (strcmp(arg, "bc-read") == 0)                      opts->bc_read = atoi(optarg);
                    else if (strcmp(arg, "sec-bc-read") == 0)                  opts->sec_bc_read = atoi(optarg);
                    else if (strcmp(arg, "first-cycle") == 0)                  ia_push(opts->first_cycle,atoi(optarg));
                    else if (strcmp(arg, "final-cycle") == 0)                  ia_push(opts->final_cycle,atoi(optarg));
                    else if (strcmp(arg, "first-index-cycle") == 0)            ia_push(opts->first_index_cycle,atoi(optarg));
                    else if (strcmp(arg, "final-index-cycle") == 0)            ia_push(opts->final_index_cycle,atoi(optarg));
                    else if (strcmp(arg, "add-cluster-index-tag") == 0)        opts->add_cluster_index_tag = true;
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

    //if (argc > 0) opts->input_name = strdup(argv[0]);
    optind = 0;

    // some validation and tidying
    if (!opts->intensity_dir) {
        fprintf(stderr,"You must specify an intensity directory (-i or --intensity-dir)\n");
        usage(stderr); return NULL;
    }

    if (opts->lane <= 0) {
        fprintf(stderr,"You must specify a lane number (-l or --lane)\n");
        usage(stderr); return NULL;
    }

    if (opts->lane > 999) {
        fprintf(stderr,"I can't handle a lane number greater than 999\n");
        usage(stderr); return NULL;
    }

    if (!opts->output_file) {
        fprintf(stderr,"You must specify an output file (-o or --output-file)\n");
        usage(stderr); return NULL;
    }

    if (opts->compression_level && !isdigit(opts->compression_level)) {
        fprintf(stderr, "compression-level must be a digit in the range [0..9], not '%c'\n", opts->compression_level);
        usage(stderr); return NULL;
    }

    // Set defaults
    if (!opts->read_group_id) opts->read_group_id = strdup("1");
    if (!opts->library_name) opts->library_name = strdup("unknown");
    if (!opts->sample_alias) opts->sample_alias = strdup(opts->library_name);
    if (!opts->sequencing_centre) opts->sequencing_centre = strdup("SC");
    if (!opts->barcode_tag) opts->barcode_tag = strdup(DEFAULT_BARCODE_TAG);
    if (!opts->quality_tag) opts->quality_tag = strdup(DEFAULT_QUALITY_TAG);
    if (!opts->bc_read) opts->bc_read = 1;
    if (!opts->sec_bc_read) opts->sec_bc_read = opts->bc_read;
    if (!opts->platform) opts->platform = strdup("ILLUMINA");

    if (!opts->run_folder) {
        // default is two levels up from intensity dierectory
        opts->run_folder = calloc(1, strlen(opts->intensity_dir) + 7);
        sprintf(opts->run_folder, "%s/../..", opts->intensity_dir);
    }

    if (!opts->basecalls_dir) {
        opts->basecalls_dir = calloc(1, strlen(opts->intensity_dir) + strlen("/BaseCalls") + 1);
        sprintf(opts->basecalls_dir, "%s/%s", opts->intensity_dir, "BaseCalls");
    }

    // rationalise directories
    char *tmp;
    tmp = opts->intensity_dir; 
    opts->intensity_dir = realpath(tmp, NULL); 
    if (!opts->intensity_dir) { perror("intensity-dir"); return NULL; }
    free(tmp);

    tmp = opts->basecalls_dir; 
    opts->basecalls_dir = realpath(tmp, NULL); 
    if (!opts->basecalls_dir) { perror("basecalls-dir"); return NULL; }
    free(tmp);

    tmp = opts->run_folder; 
    opts->run_folder = realpath(tmp, NULL); 
    if (!opts->run_folder) { perror("run_folder"); return NULL; }
    free(tmp);

    if (!opts->platform_unit) {
        // default is runfolder + lane
        char *rf = basename(opts->run_folder);
        opts->platform_unit = calloc(1, strlen(rf) + 5);
        sprintf(opts->platform_unit, "%s_%d", rf, opts->lane);
    }

    // read XML files
    opts->intensityConfig = loadXML(opts->intensity_dir, "config.xml", opts->verbose);
    opts->basecallsConfig = loadXML(opts->basecalls_dir, "config.xml", opts->verbose);
    opts->parametersConfig = loadXML(opts->run_folder, "runParameters.xml", opts->verbose);
    opts->runinfoConfig = loadXML(opts->run_folder, "RunInfo.xml", opts->verbose);


    //dumpOpts(opts);

    return opts;
}

void sam_hdr_unparse(SAM_hdr *sh, bam_hdr_t *h)
{
    free(h->text);
    sam_hdr_rebuild(sh);
    h->text = strdup(sam_hdr_str(sh));
    h->l_text = sam_hdr_length(sh);
    sam_hdr_free(sh);
}

int addHeader(samFile *output_file, bam_hdr_t *output_header, opts_t *opts)
{
    SAM_hdr *sh = sam_hdr_parse_(output_header->text,output_header->l_text);

    // Add header line
    sam_hdr_add(sh, "HD", "VN", "1.5", "SO", "unsorted", NULL, NULL);

    // Add RG line
    sam_hdr_add(sh, "RG", "ID", opts->read_group_id, 
                    "DT", opts->run_start_date,
                    "PU", opts->platform_unit,
                    "LB", opts->library_name,
                    "PG", "SCS",
                    "SM", opts->sample_alias,
                    "CN", opts->sequencing_centre,
                    "PL", opts->platform,
                    opts->study_name ? "DS": NULL, opts->study_name ? opts->study_name: NULL,
                    NULL, NULL);

    // Add PG lines
    sam_hdr_add(sh, "PG",
                    "ID", "SCS",
                    "VN", getXMLAttr(opts->intensityConfig, "/ImageAnalysis/Run/Software", "Version"), // TODO we should look first in the parameterConfig file. Ditto below
                    "PN", getXMLAttr(opts->intensityConfig, "/ImageAnalysis/Run/Software", "Name"),
                    "DS", "Controlling software on instrument",
                    NULL, NULL);

    sam_hdr_add(sh, "PG",
                    "ID", "basecalling",
                    "PP", "SCS",
                    "VN", getXMLAttr(opts->basecallsConfig, "/BaseCallAnalysis/Run/Software", "Version"),
                    "PN", getXMLAttr(opts->basecallsConfig, "/BaseCallAnalysis/Run/Software", "Name"),
                    "DS", "Basecalling Package",
                    NULL, NULL);

    sam_hdr_add(sh, "PG",
                    "ID", "bambi",
                    "PP", "basecalling",
                    "VN", bambi_version(),
                    "CL", opts->argv_list,
                    "PN", "bambi",
                    "DS", "Convert Illumina BCL to BAM or SAM file",
                    NULL, NULL);

    sam_hdr_unparse(sh,output_header);
    if (sam_hdr_write(output_file, output_header) != 0) {
        fprintf(stderr, "Could not write output file header\n");
        return 1;
    }
    return 0;
}

/*
 * load tile list from basecallsConfig or intensityConfig
 */
ia_t *getTileList(opts_t *opts)
{
    ia_t *tiles = ia_init(100);
    xmlXPathObjectPtr ptr;
    char *xpath = calloc(1,64);
    xmlDocPtr doc;

    doc = opts->basecallsConfig ? opts->basecallsConfig : opts->intensityConfig;

    sprintf(xpath, "//TileSelection/Lane[@Index=\"%d\"]/Tile", opts->lane);
    assert(strlen(xpath) < 64);
    ptr = getnodeset(doc, xpath);

    if (ptr && ptr->nodesetval) {
        int n;
        for (n=0; n < ptr->nodesetval->nodeNr; n++) {
            char * t = (char *)ptr->nodesetval->nodeTab[n]->children->content;
            if (t) ia_push(tiles,atoi(t));
        }
        xmlXPathFreeObject(ptr);
    }

    // TODO add tile range to tile array

    // TODO if (!tiles) calcTileList();

    // TODO filter tile list by command line options

    ia_sort(tiles);

    return tiles;
}

va_t *getCycleRange(opts_t *opts)
{
    va_t *cycleRange = va_init(100, freeCycleRange);
    xmlDocPtr doc;
    xmlXPathObjectPtr ptr;

    // TODO try reading from runInfo, then paramater config

    // TODO what if there is a barCodeCycleList ?

    if (va_isEmpty(cycleRange)) {
        doc = opts->basecallsConfig ? opts->basecallsConfig : opts->intensityConfig;
        ptr = getnodeset(doc, "//RunParameters/Reads");        
        if (ptr && ptr->nodesetval) {
            int n;
            for (n=0; n < ptr->nodesetval->nodeNr; n++) {
                xmlNodePtr np = ptr->nodesetval->nodeTab[n];
                char name[64];
                cycleRangeEntry_t *cr = calloc(1,sizeof(cycleRangeEntry_t));
                int readIndex = atoi((char *)xmlGetProp(np,(xmlChar *)"Index"));
                sprintf(name,"read%d",readIndex);
                cr->readname = strdup(name);
                np = np->children;
                while (np->next) {
                    if (strcmp((char*)np->name,"FirstCycle") == 0) {
                        cr->first = atoi((char*)np->children->content);
                    }
                    if (strcmp((char*)np->name,"LastCycle") == 0) {
                        cr->last = atoi((char*)np->children->content);
                    }
                    np = np->next;
                }
                va_push(cycleRange, cr);
            }
        }
    }
    return cycleRange;
}

/*
 * Open the position file
 *
 * Try looking for _pos.txt, .clocs, locs files in that order
 *
 * Open and return the first one found, or NULL if not found.
 */

FILE *openPositionFile(int tile, opts_t *opts)
{
    FILE *posFile = NULL;

    char *fname = calloc(1, strlen(opts->intensity_dir)+64);

    sprintf(fname, "%s/s_%d_%04d_pos.txt", opts->intensity_dir, opts->lane, tile);
    posFile = fopen(fname,"rb");
    if (opts->verbose && posFile) printf("Opened %s\n", fname);

    if (!posFile) {
        sprintf(fname, "%s/L%03d/s_%d_%04d.clocs", opts->intensity_dir, opts->lane, opts->lane, tile);
        posFile = fopen(fname,"rb");
        if (opts->verbose && posFile) printf("Opened %s\n", fname);
    }

    if (!posFile) {
        sprintf(fname, "%s/L%03d/s_%d_%04d.locs", opts->intensity_dir, opts->lane, opts->lane, tile);
        posFile = fopen(fname,"rb");
        if (opts->verbose && posFile) printf("Opened %s\n", fname);
    }

    free(fname);
    return posFile;

}

int processTile(int tile, samFile *output_file, bam_hdr_t *output_header, va_t *cycleRange, opts_t *opts)
{
    if (opts->verbose) printf("Processing Tile %d\n", tile);
    FILE *posFile = openPositionFile(tile, opts);
    if (!posFile) {
        fprintf(stderr,"Can't find position file for Tile %d\n", tile);
        return 1;
    }

    return 0;
}

void createBAM(samFile *output_file, bam_hdr_t *output_header, opts_t *opts)
{
    ia_t *tiles = getTileList(opts);
    va_t *cycleRange = getCycleRange(opts);
    int n;

    for (n=0; n < cycleRange->end; n++) {
        cycleRangeEntry_t *cr = (cycleRangeEntry_t *)cycleRange->entries[n];
        printf("CycleRange: %s\t%d\t%d\n", cr->readname, cr->first, cr->last);
    }

    for (n=0; n < tiles->end; n++) {
        if (processTile(tiles->entries[n], output_file, output_header, cycleRange, opts)) {
            fprintf(stderr,"Error processing tile %d\n", tiles->entries[n]);
            break;
        }
    }

    ia_free(tiles);
}

/*
 * Main code
 */
static int i2b(opts_t* opts)
{
    int retcode = 1;
    samFile *output_file = NULL;
    bam_hdr_t *output_header = NULL;
    htsFormat *out_fmt = NULL;
    char mode[] = "wbC";

    while (1) {

        /*
         * Open output file and header
         */
        if (opts->output_fmt) {
            out_fmt = calloc(1,sizeof(htsFormat));
            if (hts_parse_format(out_fmt, opts->output_fmt) < 0) {
                fprintf(stderr,"Unknown output format: %s\n", opts->output_fmt);
                break;
            }
        }
        mode[2] = opts->compression_level ? opts->compression_level : '\0';
        output_file = hts_open_format(opts->output_file, mode, out_fmt);
        free(out_fmt);
        if (!output_file) {
            fprintf(stderr, "Could not open output file (%s)\n", opts->output_file);
            break;
        }

        output_header = bam_hdr_init();
        output_header->text = calloc(1,1); output_header->l_text=0;

        if (!output_header) {
            fprintf(stderr, "Failed to initialise output header\n");
            break;
        }

        addHeader(output_file, output_header, opts);
        createBAM(output_file, output_header, opts);

        retcode = 0;
        break;
    }

    // tidy up after us
    if (output_header) bam_hdr_destroy(output_header);
    if (output_file) sam_close(output_file);
    
    return retcode;
}

/*
 * called from bambi to perform Illumina to BAM conversion
 *
 * Parse the command line arguments, then call the main i2b function
 *
 * returns 0 on success, 1 if there was a problem
 */
int main_i2b(int argc, char *argv[])
{
    int ret = 1;
    opts_t* opts = i2b_parse_args(argc, argv);
    if (opts) {
        ret = i2b(opts);
    }
    free_opts(opts);
    return ret;
}
