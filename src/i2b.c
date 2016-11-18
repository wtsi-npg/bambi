/*  i2b.c -- index i2br subcommand.

    Copyright (C) 2016 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
#include <time.h>
#include <fcntl.h>

#include <cram/sam_header.h>

#include "posfile.h"
#include "filterfile.h"
#include "bclfile.h"
#include "array.h"
#include "parse.h"

#define DEFAULT_BARCODE_TAG "BC"
#define DEFAULT_QUALITY_TAG "QT"

char *strptime(const char *s, const char *format, struct tm *tm);


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
 * BCL Read and File arrays
 */
typedef struct {
    char *readname;
    va_t *bclFileArray;
    va_t *sclFileArray;
} bclReadArrayEntry_t;

typedef struct {
    bclfile_t *bcl;
} bclFileArrayEntry_t;

void freeBCLFileArray(void *ent)
{
    bclfile_t *bcl = (bclfile_t *)ent;
    bclfile_close(bcl);
}

void freeBCLReadArray(void *ent)
{
    bclReadArrayEntry_t *ra = (bclReadArrayEntry_t *)ent;
    free(ra->readname);
    va_free(ra->bclFileArray);
    va_free(ra->sclFileArray);
    free(ra);
}

typedef struct {
    int tile;
    int clusters;
} tileIndexEntry_t;

void freetileIndexArray(void *ent)
{
    free(ent);
}

void freeTagArray(void *ent)
{
    free(ent);
}

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
    va_t *barcode_tag;
    va_t *quality_tag;
    ia_t *bc_read;
    ia_t *first_cycle;
    ia_t *final_cycle;
    ia_t *first_index_cycle;
    ia_t *final_index_cycle;
    xmlDocPtr intensityConfig;
    xmlDocPtr basecallsConfig;
    xmlDocPtr parametersConfig;
    xmlDocPtr runinfoConfig;
} opts_t;


/*
 * Release all the options
 */

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
    va_free(opts->barcode_tag);
    va_free(opts->quality_tag);
    ia_free(opts->bc_read);
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

/*
 * do something clever with an XML document
 */
xmlXPathObjectPtr getnodeset(xmlDocPtr doc, char *xpath)
{
    xmlXPathContextPtr context;
    xmlXPathObjectPtr result;

    context = xmlXPathNewContext(doc);
    if (context == NULL) {
        fprintf(stderr,"Error in xmlXPathNewContext\n");
        return NULL;
    }
    result = xmlXPathEvalExpression((xmlChar *)xpath, context);
    xmlXPathFreeContext(context);
    if (result == NULL) {
        fprintf(stderr,"Error in xmlXPathEvalExpression\n");
        return NULL;
    }
    if(xmlXPathNodeSetIsEmpty(result->nodesetval)){
        xmlXPathFreeObject(result);
        return NULL;
    }
    return result;
}
/*
 * Read an attribute for an xpath from an XML doc
 */
static char *getXMLAttr(xmlDocPtr doc, char *node, char *attr)
{
    char *v = NULL;
    xmlNodeSetPtr nodeset;
    xmlXPathObjectPtr result;

    if (!doc) return v;
    result = getnodeset(doc, node);
    if (result) {
        nodeset = result->nodesetval;
        v = (char *)xmlGetProp(nodeset->nodeTab[0], (xmlChar *)attr);
        xmlXPathFreeObject (result);
    }
    return v;
}

static int getXMLAttr_int(xmlDocPtr doc, char *node, char *attr)
{
    int n = 0;
    char *v = getXMLAttr(doc,node,attr);
    if (v) {
        n = atoi(v);
        free(v);
    }
    return n;
}

static int xmlGetProp_int(xmlNodePtr node, char *tag)
{
    int n=0;
    char *v = (char*)xmlGetProp(node,(xmlChar*)tag);
    if (v) {
        n = atoi(v);
        free(v);
    }
    return n;
}

/*
 * Read the value for an xpath for a given XML doc
 */
static char *getXMLVal(xmlDocPtr doc, char *xpath)
{
    char *val = NULL;

    if (!doc) return val;
    xmlXPathObjectPtr ptr = getnodeset(doc, xpath);

    if (ptr && ptr->nodesetval) {
        val = strdup((char *)ptr->nodesetval->nodeTab[0]->children->content);
    }
    xmlXPathFreeObject(ptr);
    return val;
}

/*
 * Load an XML file into a xmlDoc pointer
 */
xmlDocPtr loadXML(char *dir, char *fname, int verbose)
{
    xmlDocPtr doc;
    char *tmp = calloc(1, strlen(dir) + strlen(fname) + 2);
    sprintf(tmp, "%s/%s", dir, fname);
    doc = xmlReadFile(tmp, NULL, XML_PARSE_NOWARNING);
    if (!doc) {
        if (verbose) fprintf(stderr, "WARNING: Failed to parse %s/%s\n", dir, fname);
    } else {
        if (verbose) fprintf(stderr, "Opened XML file: %s/%s\n", dir, fname);
    }
    free(tmp);
    return doc;
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
"       --barcode-tag                   comma separated list of tag names for barcode sequences. [default: " DEFAULT_BARCODE_TAG "]\n"
"       --quality-tag                   comma separated list of tag name for barcode qualities. [default: " DEFAULT_QUALITY_TAG "]\n"
"       --sec-barcode-tag               DEPRECATED: Tag name for second barcode sequence. [default: null]\n"
"       --sec-quality-tag               DEPRACATED: Tag name for second barcode quality. [default: null]\n"
"       --bc-read                       comma separated list of Which reads (1 or 2) should the barcode sequences and qualities be added to?\n"
"                                       [default: 1]\n"
"       --sec-bc-read                   DEPRACATED: Which read (1 or 2) should the second barcode sequence and quality be added to?\n"
"                                       [default: bc-read]\n"
"       --first-cycle                   First cycle for each standard (non-index) read. Comma separated list.\n"
"       --final-cycle                   Last cycle for each standard (non-index) read. Comma separated list.\n"
"       --first-index-cycle             First cycle for each index read. Comma separated list.\n"
"       --final-index-cycle             Last cycle for each index read. Comma separated list.\n"
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
        { NULL, 0, NULL, 0 }
    };

    opts_t* opts = calloc(sizeof(opts_t), 1);
    if (!opts) { perror("cannot allocate option parsing memory"); return NULL; }

    opts->argv_list = stringify_argv(argc+1, argv-1);
    if (opts->argv_list[strlen(opts->argv_list)-1] == ' ') opts->argv_list[strlen(opts->argv_list)-1] = 0;

    opts->bc_read = ia_init(5);
    opts->first_cycle = ia_init(5);
    opts->final_cycle = ia_init(5);
    opts->first_index_cycle = ia_init(5);
    opts->final_index_cycle = ia_init(5);
    opts->barcode_tag = va_init(5, freeTagArray);
    opts->quality_tag = va_init(5, freeTagArray);

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
                    else if (strcmp(arg, "barcode-tag") == 0)                  parse_tags(opts->barcode_tag,optarg);
                    else if (strcmp(arg, "quality-tag") == 0)                  parse_tags(opts->quality_tag,optarg);
                    else if (strcmp(arg, "sec-barcode-tag") == 0)              parse_tags(opts->barcode_tag,optarg);
                    else if (strcmp(arg, "sec-quality-tag") == 0)              parse_tags(opts->quality_tag,optarg);
                    else if (strcmp(arg, "bc-read") == 0)                      parse_int(opts->bc_read,optarg);
                    else if (strcmp(arg, "sec-bc-read") == 0)                  parse_int(opts->bc_read,optarg);
                    else if (strcmp(arg, "first-cycle") == 0)                  parse_int(opts->first_cycle,optarg);
                    else if (strcmp(arg, "final-cycle") == 0)                  parse_int(opts->final_cycle,optarg);
                    else if (strcmp(arg, "first-index-cycle") == 0)            parse_int(opts->first_index_cycle,optarg);
                    else if (strcmp(arg, "final-index-cycle") == 0)            parse_int(opts->final_index_cycle,optarg);
                    else {
                        fprintf(stderr,"\nUnknown option: %s\n\n", arg); 
                        usage(stdout); free_opts(opts);
                        return NULL;
                    }
                    break;
        default:    fprintf(stderr,"Unknown option: '%c'\n", opt);
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
    if (va_isEmpty(opts->barcode_tag)) va_push(opts->barcode_tag,strdup(DEFAULT_BARCODE_TAG));
    if (va_isEmpty(opts->quality_tag)) va_push(opts->quality_tag,strdup(DEFAULT_QUALITY_TAG));
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
    if (!opts->intensity_dir) { 
        fprintf(stderr,"Can't open directory: %s\n", tmp);
        perror("intensity-dir"); 
        return NULL; 
    }
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
    if (!opts->parametersConfig) opts->parametersConfig = loadXML(opts->run_folder, "RunParameters.xml", opts->verbose);
    opts->runinfoConfig = loadXML(opts->run_folder, "RunInfo.xml", opts->verbose);

    if (!opts->run_start_date) {
        opts->run_start_date = getXMLVal(opts->intensityConfig, "//RunParameters/RunFolderDate");
    }
    if (!opts->run_start_date) {
        opts->run_start_date = getXMLVal(opts->parametersConfig, "//Setup/RunStartDate");
    }
    if (!opts->run_start_date) {
        opts->run_start_date = getXMLVal(opts->parametersConfig, "//RunParameters/RunStartDate");
    }

    if (!opts->run_start_date) {
        fprintf(stderr, "No run-start-date given, and none found in config files\n");
        return NULL;
    }

    // reformat date from yymmdd to YYYY-mm-dd
    if (strlen(opts->run_start_date) == 6) {
        char *tmp = calloc(1,64);
        struct tm tm;
        memset(&tm, 0, sizeof(struct tm));
        strptime(opts->run_start_date, "%y%m%d", &tm);
        strftime(tmp, 63, "%Y-%m-%dT00:00:00+0000", &tm);
        free(opts->run_start_date);
        opts->run_start_date = tmp;
    }

    return opts;
}

/*
 * convert SAM_hdr to bam_hdr
 */
void sam_hdr_unparse(SAM_hdr *sh, bam_hdr_t *h)
{
    free(h->text);
    sam_hdr_rebuild(sh);
    h->text = strdup(sam_hdr_str(sh));
    h->l_text = sam_hdr_length(sh);
    sam_hdr_free(sh);
}

/*
 * Add the header lines to the BAM file
 */
int addHeader(samFile *output_file, bam_hdr_t *output_header, opts_t *opts)
{
    SAM_hdr *sh = sam_hdr_parse_(output_header->text,output_header->l_text);
    char *version = NULL;
    char *pname = NULL;

    pname = getXMLAttr(opts->parametersConfig, "/ImageAnalysis/Run/Software", "Name");
    if (!pname) pname = getXMLAttr(opts->intensityConfig, "/ImageAnalysis/Run/Software", "Name");
    if (!pname) pname = getXMLVal(opts->parametersConfig, "//ApplicationName");
    if (!pname) { fprintf(stderr,"Can't find program name anywhere\n"); return 1; }

    version = getXMLAttr(opts->parametersConfig, "/ImageAnalysis/Run/Software", "Version");
    if (!version) version = getXMLAttr(opts->intensityConfig, "/ImageAnalysis/Run/Software", "Version");
    if (!version) version = getXMLVal(opts->parametersConfig, "//ApplicationVersion");
    if (!version) { fprintf(stderr,"Can't find program version anywhere\n"); return 1; }

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
                    (opts->study_name ? "DS": NULL), (opts->study_name ? opts->study_name: NULL),
                    NULL, NULL);

    // Add PG lines
    version = getXMLAttr(opts->parametersConfig, "/ImageAnalysis/Run/Software", "Version");
    if (!version) version = getXMLAttr(opts->intensityConfig, "/ImageAnalysis/Run/Software", "Version");
    if (!version) version = getXMLVal(opts->parametersConfig, "//Setup/ApplicationVersion");
    if (!version) { fprintf(stderr, "Can't find program version\n"); exit(1); }
    pname = getXMLAttr(opts->parametersConfig, "/ImageAnalysis/Run/Software", "Name");
    if (!pname) pname = getXMLAttr(opts->intensityConfig, "/ImageAnalysis/Run/Software", "Name");
    if (!pname) pname = getXMLVal(opts->parametersConfig, "//Setup/ApplicationName");
    if (!pname) { fprintf(stderr, "Can't find program name\n"); exit(1); }
    sam_hdr_add(sh, "PG",
                    "ID", "SCS",
                    "VN", version,
                    "PN", pname,
                    "DS", "Controlling software on instrument",
                    NULL, NULL);
    free(pname); free(version);

    version = getXMLAttr(opts->basecallsConfig, "/BaseCallAnalysis/Run/Software", "Version");
    pname = getXMLAttr(opts->basecallsConfig, "/BaseCallAnalysis/Run/Software", "Name");
    sam_hdr_add(sh, "PG",
                    "ID", "basecalling",
                    "PP", "SCS",
                    "VN", version ? version : "Unknown",
                    "PN", pname ? pname : "Unknown",
                    "DS", "Basecalling Package",
                    NULL, NULL);
    free(pname); free(version);

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
 * Get an ID - depending on what config files we have and what is in them
 *             this is either instrument_runid or computer_experiment
 */
char *getId(opts_t *opts)
{
    xmlDocPtr doc;
    char *runid = NULL;
    char *instrument = NULL;
    char *experiment = NULL;
    char *computer = NULL;
    char *id = NULL;

    doc = opts->basecallsConfig ? opts->basecallsConfig : opts->intensityConfig;

    runid = getXMLVal(doc, "//RunParameters/RunFolderId");
    instrument = getXMLVal(doc, "//RunParameters/Instrument");

    if (instrument && runid) {
        id = calloc(1, strlen(instrument) + strlen(runid) + 2);
        sprintf(id, "%s_%s", instrument, runid);
    }

    if (!id) {
        experiment = getXMLVal(opts->parametersConfig, "//Setup/ExperimentName");
        computer = getXMLVal(opts->parametersConfig, "//Setup/ComputerName");
        if (experiment && computer) {
            id = calloc(1, strlen(experiment) + strlen(computer) + 2);
            sprintf(id, "%s_%s", computer, experiment);
        }
    }

    free(instrument); free(runid); free(experiment); free(computer);

    return id;
}

/*
 * Load the tile index array (the BCI file)
 * This is only for NextSeq 
 */
va_t *getTileIndex(opts_t *opts)
{
    va_t *tileIndex = NULL;
    char *fname = calloc(1,strlen(opts->basecalls_dir)+64);
    sprintf(fname, "%s/L%03d/s_%d.bci", opts->basecalls_dir, opts->lane, opts->lane);
    int fhandle = open(fname,O_RDONLY);
    if (fhandle < 0) {
        if (opts->verbose) fprintf(stderr,"Can't open BCI file %s\n", fname);
    } else {
        tileIndex = va_init(100,freetileIndexArray);
        int n;
        do {
            tileIndexEntry_t *ti = calloc(1, sizeof(tileIndexEntry_t));
            n = read(fhandle, &ti->tile, 4);
            n = read(fhandle, &ti->clusters, 4);
            if (n == 4) {
                va_push(tileIndex,ti);
            }
        } while (n == 4);
        close(fhandle);
    }
    return tileIndex;
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
    free(xpath);

    if (ptr && ptr->nodesetval) {
        int n;
        for (n=0; n < ptr->nodesetval->nodeNr; n++) {
            char *t = (char *)ptr->nodesetval->nodeTab[n]->children->content;
            if (t) ia_push(tiles,atoi(t));
        }
        xmlXPathFreeObject(ptr);
    } else {
        // Maybe this is a NewSeq run?
        ptr = getnodeset(opts->parametersConfig, "//SelectedTiles/Tile");
        if (ptr && ptr->nodesetval) {
            int n;
            for (n=0; n < ptr->nodesetval->nodeNr; n++) {
                char *t = (char *)ptr->nodesetval->nodeTab[n]->children->content;
                char *saveptr;
                char *lane = strtok_r(t, "_", &saveptr);
                char *tileno = strtok_r(NULL, "_", &saveptr);
                if (lane && tileno) {
                    if (atoi(lane) == opts->lane) {
                        ia_push(tiles,atoi(tileno));
                    }
                }
            }
        }
    }

    //
    // If we can't find a list of tiles anywhere, try calculating them from the FlowcellLayout
    //
    if (ia_isEmpty(tiles)) {
        int numSurfaces = getXMLAttr_int(opts->runinfoConfig, "//FlowcellLayout", "SurfaceCount");
        int numSwaths = getXMLAttr_int(opts->runinfoConfig, "//FlowcellLayout", "SwathCount");
        int numTilesPerSwath = getXMLAttr_int(opts->runinfoConfig, "//FlowcellLayout", "TileCount");
        int numSectionsPerLane = getXMLAttr_int(opts->runinfoConfig, "//FlowcellLayout", "SectionPerLane");

        char *TileNamingConvention = getXMLAttr(opts->runinfoConfig, "//FlowcellLayout/TileSet", "TileNamingConvention");
        if (TileNamingConvention && strcmp(TileNamingConvention,"FiveDigit") == 0) {
            // probably a nextSeq with 5 digit tile numbers...
            if (numSurfaces && numSwaths && numTilesPerSwath && numSectionsPerLane) {
                int ispl, isur, isw, itile;
                for (ispl = 1; ispl <= numSectionsPerLane; ispl++) {
                    for (isur = 1; isur <= numSurfaces; isur++) {
                        for (isw = 1; isw <= numSwaths; isw++) {
                            for (itile = 1; itile <= numTilesPerSwath; itile++) {
                                ia_push(tiles, 10000 * isur + 1000 * ispl + 100 * isw + itile);
                            }
                        }
                    }
                }
            }
        } else {
            // 'normal' four digit tile numbers
            if (numSurfaces && numSwaths && numTilesPerSwath) {
                int isur, isw, itile;
                for (isur = 1; isur <= numSurfaces; isur++) {
                    for (isw = 1; isw <= numSwaths; isw++) {
                        for (itile = 1; itile <= numTilesPerSwath; itile++) {
                            ia_push(tiles, 1000 * isur + 100 * isw + itile);
                        }
                    }
                }
            }
        }
        free(TileNamingConvention);
    }


    if (ia_isEmpty(tiles)) return tiles;

    ia_sort(tiles);

    // Filter tile list by command line options (mainly used for testing)
    if (opts->tile_limit && opts->first_tile==0) opts->first_tile = tiles->entries[0];

    if (opts->first_tile != 0) {
        int n;
        ia_t *new_tiles = ia_init(100);
        for (n=0; n < tiles->end; n++) {
            if (tiles->entries[n] == opts->first_tile) {
                int i, tl;
                tl = opts->tile_limit ? opts->tile_limit : tiles->end;
                for (i=n; i < n+tl; i++) {
                    if (i < tiles->end) {
                        ia_push(new_tiles,tiles->entries[i]);
                    }
                }
            }
        }
        ia_free(tiles);
        tiles = new_tiles;
        if (ia_isEmpty(tiles)) {
            fprintf(stderr,"No tiles to process\n");
            exit(1);
        }
    }

    return tiles;
}

char *getCycleName(int readCount, bool isIndex)
{
    // implements naming convention used by earlier versions of Lane.java
    char *cycleName = calloc(1,16);
;
    if (isIndex) {
        if (readCount==1) { strcpy(cycleName,"readIndex"); }
        else { sprintf(cycleName,"readIndex%d",readCount); }
    } else {
        sprintf(cycleName,"read%d",readCount);
    }
    return cycleName;
}

void getCycleRangeFromFile(va_t *cycleRange, xmlDocPtr doc)
{
    xmlXPathObjectPtr ptr;
    int readCount = 1;
    int cycleCount = 1;
    int indexCount = 1;
    int n;

    if (!doc) return;
    ptr = getnodeset(doc,"/RunInfo/Run/Reads/Read");
    if (!ptr || !ptr->nodesetval) ptr = getnodeset(doc, "/RunParameters/Setup/Reads/Read");
    if (!ptr || !ptr->nodesetval) ptr = getnodeset(doc, "/RunParameters/Reads/RunInfoRead");
    if (!ptr || !ptr->nodesetval) return;   // still can't find them. Give up.

    for (n=0; n < ptr->nodesetval->nodeNr; n++) {
        xmlNodePtr np = ptr->nodesetval->nodeTab[n];
        cycleRangeEntry_t *cr = calloc(1,sizeof(cycleRangeEntry_t));
        int numCycles = xmlGetProp_int(np,"NumCycles");
        char *p = (char *)xmlGetProp(np,(xmlChar *)"IsIndexedRead");
        bool isIndexedRead = ('Y' == *p || 'y' == *p);
        free(p);
        cr->readname = getCycleName(isIndexedRead ? indexCount++ : readCount++, isIndexedRead);
        cr->first = cycleCount;
        cr->last = cycleCount + numCycles - 1;
        va_push(cycleRange,cr);
        cycleCount += numCycles;
    }
}

/*
 * Check if we have two cycle indexes, which are consecutive, and have *no* secondary barcode tag.
 * If so, merge the two indexes into one.
 */
void mergeIndexes(va_t *cycleRange, opts_t *opts)
{
    int i1=-1, i2=-1;
    int n;
    if (opts->barcode_tag->end > 1) return;     // there is a secondary barcode tag specified

    for (n = 0; n < cycleRange->end; n++) {
        cycleRangeEntry_t *cr = cycleRange->entries[n];
        if (strcmp(cr->readname, "readIndex") == 0) i1 = n;
        if (strcmp(cr->readname, "readIndex2") == 0) i2 = n;
    }

    if (i1 > -1 && i2 > -1) {
        // we have two indexes
        cycleRangeEntry_t *cr1 = cycleRange->entries[i1];
        cycleRangeEntry_t *cr2 = cycleRange->entries[i2];
        if (cr1->last+1 == cr2->first) {
            // and they look consecutive to me
            // so merge them
            cr1->last = cr2->last;
            // remove i2 entry
            cycleRange->end--;
            free(cycleRange->entries[i2]);
            for (n=i2; n < cycleRange->end; n++) {
                cycleRange->entries[n] = cycleRange->entries[n+1];
            }
        }
    }
}


/*
 * Try to find a cycle range from somewhere
 */
va_t *getCycleRange(opts_t *opts)
{
    va_t *cycleRange = va_init(100, freeCycleRange);
    xmlDocPtr doc;
    xmlXPathObjectPtr ptr = NULL;
    int n;

    //
    // read from command line options
    //
    if (!ia_isEmpty(opts->first_cycle)) {
        for (n=0; n < opts->first_cycle->end; n++) {
            cycleRangeEntry_t *cr = calloc(1,sizeof(cycleRangeEntry_t));
            cr->readname = getCycleName(n+1,false);
            cr->first = opts->first_cycle->entries[n];
            cr->last = opts->final_cycle->entries[n];
            va_push(cycleRange,cr);
        }
        for (n=0; n < opts->first_index_cycle->end; n++) {
            cycleRangeEntry_t *cr = calloc(1,sizeof(cycleRangeEntry_t));
            cr->readname = getCycleName(n+1,true);
            cr->first = opts->first_index_cycle->entries[n];
            cr->last = opts->final_index_cycle->entries[n];
            va_push(cycleRange,cr);
        }
    }

    if (va_isEmpty(cycleRange)) getCycleRangeFromFile(cycleRange, opts->runinfoConfig);
    if (va_isEmpty(cycleRange)) getCycleRangeFromFile(cycleRange, opts->parametersConfig);

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
                int readIndex = xmlGetProp_int(np,"Index");
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

    if (ptr) xmlXPathFreeObject(ptr);
    mergeIndexes(cycleRange, opts);
    return cycleRange;
}

/*
 * Find cluster number for a given tile
 * Abort if tile not found
 */
int findClusterNumber(int tile, va_t *tileIndex)
{
    int n;
    int clusterNumber = 0;
    for (n=0; n < tileIndex->end; n++) {
        tileIndexEntry_t *ti = (tileIndexEntry_t *)tileIndex->entries[n];
        if (ti->tile == tile)
            return clusterNumber;
        clusterNumber += ti->clusters;
    }
    fprintf(stderr,"findClusterNumber(%d) : no such tile\n", tile);
    exit(1);
}

int findClusters(int tile, va_t *tileIndex)
{
    int n;
    for (n=0; n < tileIndex->end; n++) {
        tileIndexEntry_t *ti = (tileIndexEntry_t *)tileIndex->entries[n];
        if (ti->tile == tile)
            return ti->clusters;
    }
    fprintf(stderr,"findClusters(%d) : no such tile\n", tile);
    exit(1);
}

/*
 * Open the position file
 *
 * Try looking for _pos.txt, .clocs, locs files in that order
 *
 * Open and return the first one found, or NULL if not found.
 */

posfile_t *openPositionFile(int tile, va_t *tileIndex, opts_t *opts)
{
    posfile_t *posfile = NULL;

    char *fname = calloc(1, strlen(opts->intensity_dir)+64);

    sprintf(fname, "%s/s_%d_%04d_pos.txt", opts->intensity_dir, opts->lane, tile);
    posfile = posfile_open(fname);
    if (opts->verbose && !posfile->errmsg) fprintf(stderr,"Opened %s\n", fname);

    if (posfile->errmsg) {
        sprintf(fname, "%s/L%03d/s_%d_%04d.clocs", opts->intensity_dir, opts->lane, opts->lane, tile);
        free(posfile->errmsg); free(posfile);
        posfile = posfile_open(fname);
        if (opts->verbose && !posfile->errmsg) fprintf(stderr,"Opened %s\n", fname);
    }

    if (posfile->errmsg) {
        sprintf(fname, "%s/L%03d/s_%d_%04d.locs", opts->intensity_dir, opts->lane, opts->lane, tile);
        free(posfile->errmsg); free(posfile);
        posfile = posfile_open(fname);
        if (opts->verbose && !posfile->errmsg) fprintf(stderr,"Opened %s\n", fname);
    }

    if (posfile->errmsg) {
        sprintf(fname, "%s/s.locs", opts->intensity_dir);
        free(posfile->errmsg); free(posfile);
        posfile = posfile_open(fname);
        if (opts->verbose && !posfile->errmsg) fprintf(stderr,"Opened %s\n", fname);
    }

    // if still not found, try NewSeq format files
    if (posfile->errmsg) {
        sprintf(fname, "%s/s_%d_pos.txt", opts->intensity_dir, opts->lane);
        posfile = posfile_open(fname);
        if (opts->verbose && !posfile->errmsg) fprintf(stderr,"Opened %s\n", fname);

        if (posfile->errmsg) {
            sprintf(fname, "%s/L%03d/s_%d.clocs", opts->intensity_dir, opts->lane, opts->lane);
            posfile = posfile_open(fname);
            if (opts->verbose && !posfile->errmsg) fprintf(stderr,"Opened %s\n", fname);
        }

        if (posfile->errmsg) {
            sprintf(fname, "%s/L%03d/s_%d.locs", opts->intensity_dir, opts->lane, opts->lane);
            posfile = posfile_open(fname);
            if (opts->verbose && !posfile->errmsg) fprintf(stderr,"Opened %s\n", fname);
        }

        if (!posfile->errmsg) {
            if (tileIndex) {
                posfile_seek(posfile,findClusterNumber(tile,tileIndex));
            } else {
                fprintf(stderr,"Trying to open %s with no tile index\n", fname);
                posfile->errmsg = strdup("Trying to open position file with no tile index");
            }
        }
    }

    free(fname);
    return posfile;

}

/*
 * find and open the filter file
 */
filter_t *openFilterFile(int tile, va_t *tileIndex, opts_t *opts)
{
    filter_t *filter = NULL;
    char *fname = calloc(1,strlen(opts->basecalls_dir)+128); // a bit arbitrary :-(

    sprintf(fname, "%s/L%03d/s_%d_%04d.filter", opts->basecalls_dir, opts->lane, opts->lane, tile);
    filter = filter_open(fname);
    if (filter->errmsg) {
        sprintf(fname, "%s/s_%d_%04d.filter", opts->basecalls_dir, opts->lane, tile);
        filter = filter_open(fname);
    }
    if (filter->errmsg) {
        sprintf(fname, "%s/L%03d/s_%d.filter", opts->basecalls_dir, opts->lane, opts->lane);
        filter = filter_open(fname);
    }

    if (opts->verbose && !filter->errmsg) fprintf(stderr,"Opened filter file %s\n", fname);

    if (tileIndex) filter_seek(filter,findClusterNumber(tile,tileIndex));

    free(fname);
    return filter;
}

/*
 * Open a single bcl (or scl) file
 */
bclfile_t *openBclFile(char *basecalls, int lane, int tile, int cycle, char *ext, va_t *tileIndex)
{
    char *fname = calloc(1, strlen(basecalls)+128);
    if (tileIndex) {    // NextSeq format
        sprintf(fname, "%s/L%03d/%04d.%s", basecalls, lane, cycle, ext);
    } else {
        sprintf(fname, "%s/L%03d/C%d.1/s_%d_%04d.%s", basecalls, lane, cycle, lane, tile, ext);
    }
    bclfile_t *bcl = bclfile_open(fname);
    if (bcl->errmsg) {
        fprintf(stderr,"Can't open %s\n%s\n", fname, bcl->errmsg);
        return NULL;
    }

    free(fname);

    if (tileIndex) bclfile_seek(bcl, findClusterNumber(tile,tileIndex));

    return bcl;
}

/*
 * Find and open all the relevant bcl and scl files
 */
va_t *openBclFiles(va_t *cycleRange, opts_t *opts, int tile, va_t *tileIndex)
{
    int n, cycle, nCycles;

    va_t *bclReadArray = va_init(5,freeBCLReadArray);

    for (n=0; n < cycleRange->end; n++) {
        cycleRangeEntry_t *cr = cycleRange->entries[n];
        bclReadArrayEntry_t *ra = calloc(1, sizeof(bclReadArrayEntry_t));
        ra->readname = strdup(cr->readname);
        nCycles = cr->last - cr->first + 1;
        ra->bclFileArray = va_init(nCycles, freeBCLFileArray);
        ra->sclFileArray = va_init(nCycles, freeBCLFileArray);

        for (cycle = cr->first; cycle <= cr->last; cycle++) {
            bclfile_t *bcl = openBclFile(opts->basecalls_dir, opts->lane, tile, cycle, "bcl", tileIndex);
            va_push(ra->bclFileArray, bcl);

            if (opts->generate_secondary_basecalls) {
                bclfile_t *bcl = openBclFile(opts->basecalls_dir, opts->lane, tile, cycle, "scl", tileIndex);
                va_push(ra->sclFileArray, bcl);
            }
        }
        va_push(bclReadArray,ra);
    }

    return bclReadArray;
}

/*
 * calculate and return the readname
 */
char *getReadName(char *id, int lane, int tile, int x, int y)
{
    char *readName = calloc(1, 128);

    if (id && *id) {
        sprintf(readName, "%s:%d:%d:%d:%d", id, lane, tile, x, y);
    } else {
        sprintf(readName, "%d:%d:%d:%d", lane, tile, x, y);
    }
    if (strlen(readName) > 127) {
        fprintf(stderr,"readName too long: %s\n", readName);
        exit(1);
    }
    return readName;
}

bool readArrayContains(va_t *bclReadArray, char *readname)
{
    int n;
    for (n=0; n < bclReadArray->end; n++) {
        bclReadArrayEntry_t *ra = bclReadArray->entries[n];
        if (strcmp(readname, ra->readname) == 0) return true;
    }
    return false;
}

/*
 * read all the bases and qualities for a given read name ("read1" or "read2")
 */
void getBases(va_t *bclReadArray, char *readname, va_t *bases, va_t *qualities, bool convert_qual)
{
    int n;
    for (n=0; n < bclReadArray->end; n++) {
        bclReadArrayEntry_t *ra = bclReadArray->entries[n];
        if (strcmp(ra->readname, readname) == 0) {
            char *b = calloc(1, ra->bclFileArray->end+1);
            char *q = calloc(1, ra->bclFileArray->end+1);
            int i;
            for (i=0; i < ra->bclFileArray->end; i++) {
                bclfile_t *bcl = ra->bclFileArray->entries[i];
                if (bclfile_next(bcl) < 0) {
                    fprintf(stderr,"Failed to read bcl file\n");
                    exit(1);
                }
                b[i] = bcl->base;
                q[i] = bcl->quality + (convert_qual ? 33 : 0);
            }
            va_push(bases,b);
            va_push(qualities,q);
            break;
        }
    }
}

/*
 * set the BAM flag
 */
int setFlag(bool second, bool filtered, bool ispaired)
{
    int flags = 0;

    flags |= BAM_FUNMAP;
    if (filtered) flags |= BAM_FQCFAIL;
    if (ispaired) {
        flags |= BAM_FPAIRED;
        flags |= BAM_FMUNMAP;
        if (second) flags |= BAM_FREAD2;
        else        flags |= BAM_FREAD1;
    }
    return flags;
}

/*
 * Write a BAM record
 */
void writeRecord(int flags, opts_t *opts, char *readName, 
                 char *bases, char *qualities, va_t *ib, va_t *iq,
                 samFile *output_file, bam_hdr_t *output_header)
{
    bam1_t *bam = bam_init1();
    int n;

    int r = bam_construct_seq(&bam, 0, readName, strlen(readName),
                                flags, -1, 0, 0, 0, 0, (uint32_t*)"", -1, 0, 0, strlen(bases), bases, qualities);
    if (r) {
        fprintf(stderr,"bam_construct_seq() failed\n");
        exit(1);
    }

    // add read group
    bam_aux_append(bam, "RG", 'Z', strlen(opts->read_group_id)+1, (uint8_t *)opts->read_group_id);

    if (ib->end > opts->barcode_tag->end) {
        fprintf(stderr, "Not enough barcode tags. This is probably a dual index run with only one barcode tag specified\n");
        exit(1);
    }

    // add index tags
    for (n=0; n < ib->end; n++) {
        char *tag = ib->entries[n];
        bam_aux_append(bam, opts->barcode_tag->entries[n], 'Z', strlen(tag)+1, (uint8_t *)tag);
        tag = iq->entries[n];
        bam_aux_append(bam, opts->quality_tag->entries[n], 'Z', strlen(tag)+1, (uint8_t *)tag);
    }

    r = sam_write1(output_file, output_header, bam);
    if (r <= 0) {
        fprintf(stderr, "Problem writing record %s  : r=%d\n", readName,r);
        exit(1);
    }
    bam_destroy1(bam);
}

/*
 * Write all the BAM records for a given tile
 */
int processTile(int tile, samFile *output_file, bam_hdr_t *output_header, va_t *cycleRange, va_t *tileIndex, opts_t *opts)
{
    va_t *bclReadArray;
    int filtered;
    int max_cluster = 0;
    int nRecords = 0;
    int c;

    if (opts->verbose) fprintf(stderr,"Processing Tile %d\n", tile);
    posfile_t *posfile = openPositionFile(tile, tileIndex, opts);
    if (posfile->errmsg) {
        fprintf(stderr,"Can't find position file for Tile %d\n%s\n", tile, posfile->errmsg);
        return 1;
    }

    filter_t *filter = openFilterFile(tile,tileIndex,opts);
    if (filter->errmsg) {
        fprintf(stderr,"Can't find filter file for tile %d\n%s\n", tile, filter->errmsg);
        return 1;
    }

    if (tileIndex) max_cluster = findClusters(tile, tileIndex);

    bclReadArray = openBclFiles(cycleRange, opts, tile, tileIndex);
    char *id = getId(opts);

    bool ispaired = readArrayContains(bclReadArray, "read2");

    //
    // write all the records
    //
    while ( (filtered = filter_next(filter)) >= 0) {
        if (tileIndex && filter->current_cluster > max_cluster) break;
        filtered = !filtered;   // don't ask
        posfile_next(posfile);
        char *readName = getReadName(id, opts->lane, tile, posfile->x, posfile->y);
        va_t *bases, *qualities, *bases_index, *qualities_index, *bases_index2, *qualities_index2;
        bases = va_init(2,free); qualities = va_init(2,free);
        bases_index = va_init(5,free); qualities_index = va_init(5,free); bases_index2 = va_init(5,free); qualities_index2 = va_init(5,free);

        getBases(bclReadArray, "read1", bases, qualities, false);
        if (ispaired) getBases(bclReadArray, "read2", bases, qualities, false);

        // read each index and put into first or second read
        for (c=0; c < cycleRange->end; c++) {
            char *cname = getCycleName(c+1,true);
            if (readArrayContains(bclReadArray,cname)) {
                if (c >= opts->bc_read->end) ia_push(opts->bc_read,1);   // supply a default
                if (opts->bc_read->entries[c] == 2) getBases(bclReadArray, cname, bases_index2, qualities_index2, true);
                else                                getBases(bclReadArray, cname, bases_index, qualities_index, true);
            }
            free(cname);
        }

        if (opts->no_filter || !filtered) {
            int flags;
            flags = setFlag(false,filtered,ispaired);
            writeRecord(flags, opts, readName, bases->entries[0], qualities->entries[0], bases_index, qualities_index, output_file, output_header);
            if (ispaired) {
                flags = setFlag(true,filtered,ispaired);
                writeRecord(flags, opts, readName, bases->entries[1], qualities->entries[1], bases_index2, qualities_index2, output_file, output_header);
            }
            nRecords++;
        }

        va_free(bases); va_free(qualities);
        va_free(bases_index); va_free(qualities_index);
        va_free(bases_index2); va_free(qualities_index2);
        free(readName);
    }

    free(id);
    va_free(bclReadArray);
    filter_close(filter);
    posfile_close(posfile);

    if (opts->verbose) fprintf(stderr,"%d records written\n", nRecords);

    return 0;
}

/*
 * process all the tiles and write all the BAM records
 */
int createBAM(samFile *output_file, bam_hdr_t *output_header, opts_t *opts)
{
    int retcode = 0;
    ia_t *tiles = getTileList(opts);
    va_t *cycleRange = getCycleRange(opts);;
    va_t *tileIndex = getTileIndex(opts);

    int n;

    for (n=0; n < cycleRange->end; n++) {
        cycleRangeEntry_t *cr = (cycleRangeEntry_t *)cycleRange->entries[n];
        if (opts->verbose) fprintf(stderr,"CycleRange: %s\t%d\t%d\n", cr->readname, cr->first, cr->last);
    }

    if (tiles->end == 0) fprintf(stderr, "There are no tiles to process\n");

    for (n=0; n < tiles->end; n++) {
        if (processTile(tiles->entries[n], output_file, output_header, cycleRange, tileIndex, opts)) {
            fprintf(stderr,"Error processing tile %d\n", tiles->entries[n]);
            retcode = 1;
            break;
        }
    }

    va_free(cycleRange);
    ia_free(tiles);
    return retcode;
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

        if (addHeader(output_file, output_header, opts) != 0) {
            fprintf(stderr,"Failed to write header\n");
            break;
        }

        retcode = createBAM(output_file, output_header, opts);
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
    if (opts) ret = i2b(opts);
    free_opts(opts);
    return ret;
}
