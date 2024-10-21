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
#include <math.h>

#include "zlib.h"

#include "wqueue.h"
#include "wvec.h"
#include "wzmisc.h"

#include "hts.h"
#include "sam.h"

#include "khashl.h"

#include "coverage.h"

// Record related structs and functions
#define RECORD_QUEUE_END -2
#define RECORD_SLOT_OBSOLETE -1

KHASHL_MAP_INIT(KH_LOCAL, covg_map, cm, uint32_t, uint32_t, kh_hash_uint32, kh_eq_generic)

typedef struct regions_t {
    char     *chrm;    /* chromosome */
    size_t    cap;     /* array capacity for both starts and widths */
    size_t    n;       /* number of regions */
    uint32_t *starts; /* region starts */
    uint32_t *widths; /* region widths */
} regions_t;

DEFINE_VECTOR(regions_v, regions_t)

void destroy_regions(regions_v *regions) {
    uint32_t i;
    for (i=0; i<regions->size; i++) {
        regions_t *reg = ref_regions_v(regions, i);
        free(reg->starts);
        free(reg->widths);
        free(reg->chrm);
    }
    free_regions_v(regions);
}

void realloc_regions(regions_t *reg) {
    // Could potentially make this a parameter that can be toggled to increase speed further
    // Larger value = fewer allocations, but more excess memory used
    // Smaller value = more allocations, but less memory overhead
    size_t additional = 10000;

    reg->starts = realloc(reg->starts, (reg->n+additional)*sizeof(uint32_t));
    if (!reg->starts) {
        fprintf(stderr, "Failed to properly allocate space for region starts\n");
        exit(1);
    }

    reg->widths = realloc(reg->widths, (reg->n+additional)*sizeof(uint32_t));
    if (!reg->widths) {
        fprintf(stderr, "Failed to properly allocate space for region widths\n");
        exit(1);
    }

    reg->cap += additional;
}

// Get all regions from one chromosome
static inline regions_t *get_regions(regions_v *regions, char *chrm) {
    uint32_t i;
    regions_t *reg;
    for (i=0; i<regions->size; ++i) {
        reg = ref_regions_v(regions, i);
        if (strcmp(reg->chrm, chrm) == 0) { return reg; }
    }

    return NULL;
}

static inline regions_t *get_n_insert_region(regions_v *regions, char *chrm) {
    regions_t *reg = get_regions(regions, chrm);
    if (!reg) {
        reg = next_ref_regions_v(regions);
        reg->chrm = strdup(chrm);
        reg->starts = calloc(1, sizeof(uint32_t));
        reg->widths = calloc(1, sizeof(uint32_t));
        reg->n = 0;
        reg->cap = 1;
    }
    return reg;
}

#define regions_test(regs, i) regs[(i)>>3]&(1<<((i)&0x7))
#define regions_set(regs, i) regs[(i)>>3] |= 1<<((i)&0x7)

// Information stored for each window
typedef struct {
    int64_t   block_id; /* ID of block processed by thread */
    covg_map  *all;     /* hash map : (key: coverage, value: number of bases with coverage) */
    covg_map  *q40;     /* hash map : (key: coverage, value: number of bases with coverage) */
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
    char             *bam_fn;      /* BAM filename */
    regions_v        *cpg_regions; /* vector of CpGs */
    wqueue_t(window) *q;           /* window queue */
    wqueue_t(record) *rq;          /* records queue */
    covg_conf_t      *conf;        /* config variables */
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
    char *outfn;
    target_v *targets;
    covg_conf_t *conf;
} writer_conf_t;

static void merge(covg_map *merge_from, covg_map *merge_into, uint32_t *numerator, uint32_t *denominator) {
    khint_t k;
    int absent;
    kh_foreach(merge_from, k) {
        uint32_t covg = kh_key(merge_from, k);
        uint32_t base = kh_val(merge_from, k);

        *numerator += covg * base;
        *denominator += base;

        khint_t key = cm_put(merge_into, covg, &absent);
        if (absent) {
            kh_val(merge_into, key) = base;
        } else {
            kh_val(merge_into, key) += base;
        }
    }
}

static void process_coverage_results(covg_map *cm, uint32_t mean_numerator, uint32_t denominator, char *covg_fname, char *covdist_tag, FILE *cv_file, char *cv_tag) {
    double mean = (double)mean_numerator / (double)denominator;

    FILE *out = fopen(covg_fname, "w");
    fprintf(out, "BISCUITqc Depth Distribution - %s\ndepth\tcount\n", covdist_tag);

    khint_t k;
    uint32_t variance_numerator = 0;
    kh_foreach(cm, k) {
        uint32_t covg = kh_key(cm, k);
        uint32_t base = kh_val(cm, k);

        variance_numerator += base * (covg - mean) * (covg - mean);

        fprintf(out, "%u\t%u\n", covg, base);
    }

    fflush(out);
    fclose(out);

    double sigma = sqrt((double)variance_numerator / (double)denominator);
    fprintf(cv_file, "%s\t%lf\t%lf\t%lf\n", cv_tag, mean, sigma, sigma/mean);
}

typedef struct {
    covg_map *all_base; /* all reads base coverage */
    covg_map *q40_base; /* q40 reads cpg coverage */
    //covg_map *all_cpg; /* all reads base coverage */
    //covg_map *q40_cpg; /* q40 reads cpg coverage */
} maps_t;

static inline maps_t *init_maps() {
    maps_t *out = calloc(1, sizeof(maps_t));

    out->all_base = cm_init();
    out->q40_base = cm_init();
    //out->all_cpg  = cm_init();
    //out->q40_cpg  = cm_init();

    return out;
}

static inline void destroy_maps(maps_t *maps) {
    //cm_destroy(maps->q40_cpg);
    //cm_destroy(maps->all_cpg);
    cm_destroy(maps->q40_base);
    cm_destroy(maps->all_base);

    free(maps);
}

static void *coverage_write_func(void *data) {
    writer_conf_t *c = (writer_conf_t*) data;

    FILE *out;
    if (c->outfn) {
        out = fopen(c->outfn, "w");
    } else {
        out = stdout;
    }

    fprintf(out, "BISCUITqc Depth Distribution - ...\ndepth\tcount\n");

    maps_t *maps = init_maps();

    uint32_t num_all = 0; // coverage * (N bases with coverage)
    uint32_t den_all = 0; // N bases with coverage
    uint32_t num_q40 = 0; // coverage * (N bases with coverage)
    uint32_t den_q40 = 0; // N bases with coverage
    while (1) {
        record_t rec;
        wqueue_get(record, c->q, &rec);
        if(rec.block_id == RECORD_QUEUE_END) break;

        merge(rec.all, maps->all_base, &num_all, &den_all);
        merge(rec.q40, maps->q40_base, &num_q40, &den_q40);

        cm_destroy(rec.all);
        cm_destroy(rec.q40);
    }

    FILE *cv_table = fopen("cv_table.txt", "w");
    fprintf(cv_table, "BISCUITqc Uniformity Table\ngroup\tmu\tsigma\tcv\n");

    process_coverage_results(maps->all_base, num_all, den_all, "covdist_all_base_table.txt", "All Bases", cv_table, "all_base");
    process_coverage_results(maps->q40_base, num_q40, den_q40, "covdist_q40_base_table.txt", "Q40 Bases", cv_table, "q40_base");

    fflush(cv_table);
    fclose(cv_table);

    destroy_maps(maps);

    if (c->outfn) {
        // For stdout, will close at the end of main
        fflush(out);
        fclose(out);
    }

    return 0;
}

static void format_coverage_data(covg_map *all, covg_map *q40, uint32_t *all_covgs, uint32_t *q40_covgs, uint32_t arr_len, uint8_t *cpgs) {
    int absent_all, absent_q40;
    khint_t k_all, k_q40;

    uint32_t i;
    for (i=0; i<arr_len; i++) {
        if (cpgs && regions_test(cpgs, i)) {
            fprintf(stderr, "cpg covered base! %u\n", i);
        }

        uint8_t is_match_all = 0;
        uint8_t is_match_q40 = 0;

        // If the current coverage is the same as the last one, we know the coverage has been seen before
        // and the correct bucket is already loaded up, so shortcircuit by automatically incrementing value
        if (i > 0) {
            is_match_all = all_covgs[i] == all_covgs[i-1];
            is_match_q40 = all_covgs[i] == all_covgs[i-1];

            if (is_match_all) {
                kh_val(all, k_all) += 1;
            }
            if (is_match_q40) {
                kh_val(q40, k_q40) += 1;
            }
        }

        if (!is_match_all) {
            k_all = cm_put(all, all_covgs[i], &absent_all);
            if (absent_all) {
                kh_val(all, k_all) = 1;
            } else {
                kh_val(all, k_all) += 1;
            }
        }

        if (!is_match_q40) {
            k_q40 = cm_put(q40, q40_covgs[i], &absent_q40);
            if (absent_q40) {
                kh_val(q40, k_q40) = 1;
            } else {
                kh_val(q40, k_q40) += 1;
            }
        }
    }
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

        char *chrm = header->target_name[w.tid];

        uint8_t *cpgs = NULL;
        if (res->cpg_regions) {
            cpgs = calloc((w.end - w.beg)/8 + 1, sizeof(uint8_t));
            regions_t *reg = get_regions(res->cpg_regions, chrm);

            // If chromosome is found in CpG regions file
            int j;
            if (reg) {
                for (j=0; j<reg->n; j++) {
                    uint32_t start = reg->starts[j];
                    uint32_t width = reg->widths[j];
                    if (start+width >= w.beg && start < w.end) {
                        int k;
                        for (k=0; k<width; k++) {
                            if (start+k >= w.beg && start+k < w.end) {
                                regions_set(cpgs, start+k-w.beg);
                            }
                        }
                    }
                }
            }
        }

        uint32_t *all_covgs = calloc(w.end - w.beg, sizeof(uint32_t));
        uint32_t *q40_covgs = calloc(w.end - w.beg, sizeof(uint32_t));

        hts_itr_t *iter = sam_itr_queryi(idx, w.tid, w.beg>1?(w.beg-1):1, w.end);
        bam1_t *b = bam_init1();
        int ret;
        while ((ret = sam_itr_next(in, iter, b))>0) {
            bam1_core_t *c = &b->core;

            // Read-based filtering
            // TODO: If I want to get the exact same results as from bedtools, I
            //       will likely need to adjust my filters to match bedtools:
            //       genomeCoverageBed/genomeCoverageBed.cpp:L303 and following
            if (c->flag > 0) { // only when any flag is set
                if (c->flag & BAM_FUNMAP) continue;
                if (c->flag & BAM_FSECONDARY) continue;
                if (c->flag & BAM_FQCFAIL) continue;
                if (c->flag & BAM_FDUP) continue;
                if (c->flag & BAM_FSUPPLEMENTARY) continue;
            }

            uint32_t rpos = c->pos + 1; // 1-based reference position
            uint32_t qpos = 0;

            int i;
            for (i=0; i<c->n_cigar; ++i) {
                uint32_t op    = bam_cigar_op(bam_get_cigar(b)[i]);
                uint32_t oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
                switch(op) {
                    case BAM_CMATCH:
                    case BAM_CEQUAL:
                    case BAM_CDIFF:
                        for (j=0; j<oplen; ++j) {
                            if (w.beg <= rpos && rpos+j < w.end) {
                                uint32_t idx = rpos + j - w.beg;
                                all_covgs[idx] += 1;
                                if (c->qual >= 40) {
                                    q40_covgs[idx] += 1;
                                }
                            }
                        }
                        rpos += oplen;
                        qpos += oplen;
                        break;
                    case BAM_CINS:
                    case BAM_CSOFT_CLIP:
                        qpos += oplen;
                        break;
                    case BAM_CDEL:
                    case BAM_CREF_SKIP:
                        rpos += oplen;
                        break;
                    case BAM_CHARD_CLIP:
                    case BAM_CPAD:
                        break;
                    default:
                        fprintf(stderr, "Unknown cigar %u\n", op);
                        abort();
                }
            }
        }

        // produce coverage output
        rec.all = cm_init();
        rec.q40 = cm_init();
        format_coverage_data(rec.all, rec.q40, all_covgs, q40_covgs, w.end-w.beg, cpgs);

        // set record block id
        rec.block_id = w.block_id;
        // put output string to output queue
        wqueue_put2(record, res->rq, rec);

        bam_destroy1(b);
        hts_itr_destroy(iter);

        free(q40_covgs);
        free(all_covgs);
        free(cpgs);
    }

    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    hts_close(in);

    return 0;
}

regions_v *bed_init_regions(char *bed_fn) {
    regions_v *regions = init_regions_v(2);
    kstring_t line;
    line.l = line.m = 0;
    line.s = 0;

    // Read BED file
    regions_t *reg = 0;
    char *tok;

    gzFile fh = gzopen(bed_fn, "r");
    if (fh == NULL) {
        free(line.s);
        fprintf(stderr, "Could not find regions BED file: %s\n", bed_fn);
        gzclose(fh);
        exit(1);
    }

    uint32_t n_lines = 0;
    uint8_t first_char = 1;
    while (1) {
        int c = gzgetc(fh);
        if (c < 0 && first_char) {
            free(line.s);
            fprintf(stderr, "Regions BED file (%s) is empty\n", bed_fn);
            gzclose(fh);
            exit(1);
        }
        first_char = 0;

        if (c == '\n' || c == EOF || c < 0) {
            if (strcount_char(line.s, '\t') == 2) {
                // Get chromosome
                tok = strtok(line.s, "\t");
                if (!reg || strcmp(reg->chrm, tok) != 0) {
                    fprintf(stderr, "new chromosome found: %s\n", tok);
                    n_lines = 0;
                    reg = get_n_insert_region(regions, tok);
                }

                // Adjust number of elements
                if (n_lines == reg->cap) {
                    realloc_regions(reg);
                }

                // Get start
                tok = strtok(NULL, "\t");
                ensure_number(tok);
                uint32_t start = (uint32_t)atoi(tok);

                tok = strtok(NULL, "\t");
                ensure_number(tok);
                uint32_t end = (uint32_t)atoi(tok);

                reg->starts[reg->n] = start;
                reg->widths[reg->n] = end - start;

                n_lines++;
                reg->n++;
            }

            line.l = 0;
            if (c == EOF || c < 0) {
                break;
            }
        } else {
            kputc(c, &line);
        }
    }

    gzclose(fh);
    free(line.s);

    return regions;
}

void covg_conf_init(covg_conf_t *conf) {
    conf->step = 100000;
    conf->n_threads = 1;
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
    char *cpg_bed_fn = 0;
    char *out_fn = 0;

    covg_conf_t conf;
    covg_conf_init(&conf);

    int c;
    if (argc < 2) { return usage(); }
    while ((c=getopt(argc, argv, ":@:m:o:r:s:")) >= 0) {
        switch (c) {
            case '@': conf.n_threads = atoi(optarg); break;
            case 'o': out_fn = optarg; break;
            case 'r': region_bed_fn = optarg; break;
            case 'm': cpg_bed_fn = optarg; break;
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
    regions_v *cpg_regions = cpg_bed_fn ? bed_init_regions(cpg_bed_fn) : NULL;
    fprintf(stderr, "Number of chromosomes: %u\n", cpg_regions->size);

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
        .targets = targets,
        .conf = &conf,
    };
    pthread_create(&writer, NULL, coverage_write_func, &writer_conf);
    for (i=0; i<conf.n_threads; ++i) {
        results[i].q = wq;
        results[i].rq = writer_conf.q;
        results[i].cpg_regions = cpg_regions;
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

    if (cpg_regions)
        destroy_regions(cpg_regions);

    return 0;
}
