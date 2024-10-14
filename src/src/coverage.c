/* find coverage metrics from a bam file
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2024 Jacob.Morrison@vai.org
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "wqueue.h"
#include "wvec.h"

#include "hts.h"
#include "sam.h"

#include "coverage.h"

// Record related structs and functions
#define RECORD_QUEUE_END -2
#define RECORD_SLOT_OBSOLETE -1

typedef struct regions_t {
    char     *chrm;    /* chromosome */
    size_t    n;       /* number of reginos */
    uint32_t *locs[2]; /* region starts (idx: 0) and lengths (idx: 1) */
} regions_t;

DEFINE_VECTOR(regions_v, regions_t)

// Information stored for each window
typedef struct {
    int64_t   block_id; /* ID of block processed by thread */
    kstring_t s;        /* stores entries to print */
    int       tid;      /* contig ID number */
} record_t;

DEFINE_VECTOR(record_v, record_t)
DEFINE_WQUEUE(record, record_t)

// Window blocks for processing regions
typedef struct {
    int64_t  block_id; /* ID number for window */
    int32_t  tid;      /* contig ID number of region */
    uint32_t beg;      /* beginning of region for window */
    uint32_t end;      /* end of region for window */
} window_t;

DEFINE_WQUEUE(window, window_t)

// Shared information across threads
typedef struct {
    char             *bam_fn;  /* BAM filename */
    regions_v        *regions; /* vector of regions */
    wqueue_t(window) *q;       /* window queue */
    wqueue_t(record) *rq;      /* records queue */
    covg_conf_t      *conf;    /* config variables */
} result_t;

// Contig info
typedef struct {
    int32_t   tid;  /* contig ID number */
    char     *name; /* contig name */
    uint32_t  len;  /* contig length */
} target_t;

DEFINE_VECTOR(target_v, target_t);

static inline int compare_targets(const void *a, const void *b) {
    return strcmp(((target_t*)a)->name, ((target_t*)b)->name);
}

typedef struct {
    wqueue_t(record) *q;
    char **bam_fns;
    char *outfn;
    char *header;
    target_v *targets;
    covg_conf_t *conf;
} writer_conf_t;

void pop_record_by_block_id(record_v *records, int64_t block_id, record_t *record) {
    uint64_t i;
    record_t *r;
    for (i=0; i<records->size; ++i) {
        r = ref_record_v(records, i);
        if (r->block_id == block_id) {
            *record = *r;             /* copy the record and set slot on shelf to OBSOLETE */
            r->block_id = RECORD_SLOT_OBSOLETE;
            return;
        }
    }
    record->block_id = RECORD_SLOT_OBSOLETE;
}

void put_into_record_v(record_v *records, record_t rec) {
    uint64_t i;
    record_t *r;

    /* fill blanks */
    for (i=0; i<records->size; ++i) {
        r = ref_record_v(records, i);
        if (r->block_id == RECORD_SLOT_OBSOLETE) {
            *r = rec;
            return;
        }
    }

    /* get a new slot */
    r = next_ref_record_v(records);
    *r = rec;
    return;
}

static void *coverage_write_func(void *data) {
    writer_conf_t *c = (writer_conf_t*) data;

    FILE *out;
    if (c->outfn) {
        out = fopen(c->outfn, "w");
    } else {
        out = stdout;
    }

    int64_t next_block = 0;
    record_v *records = init_record_v(20);

    while (1) {
        record_t rec;
        wqueue_get(record, c->q, &rec);
        if(rec.block_id == RECORD_QUEUE_END) break;

        if (rec.block_id == next_block) {
            do {
                if (rec.s.s) fputs(rec.s.s, out);
                free(rec.s.s);

                // Get next block from shelf if available else return OBSOLETE and retrieve new block from queue
                next_block++;
                pop_record_by_block_id(records, next_block, &rec);
            } while (rec.block_id != RECORD_SLOT_OBSOLETE);
        } else {
            // Shelf the block if not next
            put_into_record_v(records, rec);
        }
    }

    free_record_v(records);
    if (c->outfn) {
        // For stdout, will close at the end of main
        fflush(out);
        fclose(out);
    }

    return 0;
}

static void *process_func(void *data) {
    result_t *res  = (result_t*) data;
    covg_conf_t *conf = (covg_conf_t*) res->conf;

    htsFile   *in  = hts_open(res->bam_fn, "rb");
    hts_idx_t *idx = sam_index_load(in, res->bam_fn);
    if (!idx) {
        fprintf(stderr, "[%s:%d] BAM %s is not indexed?\n", __func__, __LINE__, res->bam_fn);
        fflush(stderr);
        exit(1);
    }
    bam_hdr_t *header = sam_hdr_read(in);

    uint32_t j;

    record_t rec;
    memset(&rec, 0, sizeof(record_t));
    window_t w;
    while (1) {
        wqueue_get(window, res->q, &w);
        if (w.tid == -1) break;

        rec.tid = w.tid;
        char *chrm = header->target_name[w.tid];

        // The output string
        rec.s.l = rec.s.m = 0; rec.s.s = 0;

        ksprintf(&rec.s, "%lli\t%i\t%u\t%u\n", w.block_id, w.tid, w.beg, w.end);

        // set record block id
        rec.block_id = w.block_id;
        // put output string to output queue
        wqueue_put2(record, res->rq, rec);
    }

    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    hts_close(in);

    return 0;
}

void covg_conf_init(covg_conf_t *conf) {
    conf->step = 100000;
    conf->n_threads = 3;
}

static int usage() {
    covg_conf_t conf;
    covg_conf_init(&conf);

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: coverage [options] <in.bam>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -o STR    Output filename\n");
    fprintf(stderr, "    -r STR    BED file of regions of interest\n");
    fprintf(stderr, "    -s INT    Step size of windows [%d]\n", conf.step);
    fprintf(stderr, "    -@ INT    Number of threads [%d]\n", conf.n_threads);
    fprintf(stderr, "\n");

    return 1;
}

int main_coverage(int argc, char *argv[]) {
    char *region_bed_fn = 0;
    char *out_fn = 0;

    covg_conf_t conf;
    covg_conf_init(&conf);

    int c;
    if (argc < 2) { return usage(); }
    while ((c=getopt(argc, argv, ":@:o:r:s:")) >= 0) {
        switch (c) {
            case '@': conf.n_threads = atoi(optarg); break;
            case 'o': out_fn = optarg; break;
            case 'r': region_bed_fn = optarg; break;
            case 's': conf.step = atoi(optarg); break;
            case ':': usage(); fprintf(stderr, "Option needs an argument: -%c\n", optopt); return 1;
            case '?': usage(); fprintf(stderr, "Unrecognized option: -%c\n", optopt); return 1;
            default: return usage();
        }
    }

    if (optind + 1 > argc) {
        usage();
        fprintf(stderr, "BAM input is missing\n");
        return 1;
    }
    char *infn = argv[optind++];

    // TODO: Placeholder until I implement region coverages
    regions_v *regions = NULL;

    htsFile *in = hts_open(infn, "rb");
    if (in == NULL) {
        fprintf(stderr, "%s unable to be opened\n", infn);
        return 1;
    }
    bam_hdr_t *header = sam_hdr_read(in);

    // sort sequence name by alphabetic order, chr1, chr10, chr11 ...
    target_v *targets = init_target_v(50);
    target_t *t;
    int i;
    for (i=0; i<header->n_targets; ++i) {
        t = next_ref_target_v(targets);
        t->tid = i;
        t->name = header->target_name[i];
        t->len = header->target_len[i];
    }

    qsort(targets->buffer, targets->size, sizeof(target_t), compare_targets);
    
    // Setup multithreading
    wqueue_t(window) *wq = wqueue_init(window, conf.step);
    pthread_t *processors = calloc(conf.n_threads, sizeof(pthread_t));
    result_t *results = calloc(conf.n_threads, sizeof(result_t));

    // Setup writer
    pthread_t writer;
    writer_conf_t writer_conf = {
        .q = wqueue_init(record, conf.step),
        .outfn = out_fn,
        .header = 0,
        .targets = targets,
        .conf = &conf,
    };
    pthread_create(&writer, NULL, coverage_write_func, &writer_conf);
    for (i=0; i<conf.n_threads; ++i) {
        results[i].q = wq;
        results[i].rq = writer_conf.q;
        results[i].bam_fn = infn;
        results[i].conf = &conf;
        pthread_create(&processors[i], NULL, process_func, &results[i]);
    }

    window_t w;
    memset(&w, 0, sizeof(window_t));
    uint32_t wbeg;
    int64_t block_id=0;

    // process bam
    uint32_t j;
    for (j=0; j<targets->size; ++j) {
        t = ref_target_v(targets, j);
        for (wbeg=1; wbeg<t->len; wbeg += conf.step, block_id++) {
            w.tid = t->tid;
            w.block_id = block_id;
            w.beg = wbeg;
            w.end = wbeg + conf.step;
            if (w.end > t->len) w.end = t->len;
            wqueue_put(window, wq, &w);
        }
    }
    for (i=0; i<conf.n_threads; ++i) {
        w.tid = -1;
        wqueue_put(window, wq, &w);
    }

    for (i=0; i<conf.n_threads; ++i) {
        pthread_join(processors[i], NULL);
    }

    record_t rec = { .block_id = RECORD_QUEUE_END };
    wqueue_put2(record, writer_conf.q, rec);
    pthread_join(writer, NULL);

    wqueue_destroy(record, writer_conf.q);
    free_target_v(targets);
    free(results);
    free(processors);
    wqueue_destroy(window, wq);
    hts_close(in);
    bam_hdr_destroy(header);

    //if (regions)
    //    destroy_regions(regions);

    return 0;
}
