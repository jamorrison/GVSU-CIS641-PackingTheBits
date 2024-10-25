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

typedef struct {
    uint32_t num; /* fraction numerator */
    uint32_t den; /* fraction denominator */
} fraction_t;

static inline fraction_t init_fraction() {
    fraction_t out;
    out.num = 0;
    out.den = 0;

    return out;
}

static inline double divide(fraction_t frac) {
    return (double)frac.num / (double)frac.den;
}

typedef struct {
    fraction_t all_base;
    fraction_t q40_base;
    fraction_t all_cpg;
    fraction_t q40_cpg;
    fraction_t all_base_top;
    fraction_t q40_base_top;
} total_coverage_t;

static inline total_coverage_t init_total_coverage() {
    total_coverage_t out;

    out.all_base     = init_fraction();
    out.q40_base     = init_fraction();
    out.all_cpg      = init_fraction();
    out.q40_cpg      = init_fraction();
    out.all_base_top = init_fraction();
    out.q40_base_top = init_fraction();

    return out;
}

typedef struct {
    covg_map *all_base;     /* all reads base coverage */
    covg_map *q40_base;     /* q40 reads cpg coverage */
    covg_map *all_cpg;      /* all reads base coverage */
    covg_map *q40_cpg;      /* q40 reads cpg coverage */
    covg_map *all_base_top; /* all reads top GC content coverage */
    covg_map *q40_base_top; /* q40 reads top GC content coverage */
} maps_t;

static inline maps_t *init_maps() {
    maps_t *out = calloc(1, sizeof(maps_t));

    out->all_base     = cm_init();
    out->q40_base     = cm_init();
    out->all_cpg      = cm_init();
    out->q40_cpg      = cm_init();
    out->all_base_top = cm_init();
    out->q40_base_top = cm_init();

    return out;
}

static inline void destroy_maps(maps_t *maps) {
    cm_destroy(maps->q40_base_top);
    cm_destroy(maps->all_base_top);
    cm_destroy(maps->q40_cpg);
    cm_destroy(maps->all_cpg);
    cm_destroy(maps->q40_base);
    cm_destroy(maps->all_base);

    free(maps);
}

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
    int64_t  block_id; /* ID of block processed by thread */
    maps_t  *maps;     /* coverage hash maps */
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
    regions_v        *top_regions; /* vector of top GC content regions */
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
    uint8_t has_top;
    wqueue_t(record) *q;
    char *prefix;
    target_v *targets;
    covg_conf_t *conf;
} writer_conf_t;

static void merge(covg_map *merge_from, covg_map *merge_into, fraction_t *frac) {
    khint_t k;
    int absent;
    kh_foreach(merge_from, k) {
        uint32_t covg = kh_key(merge_from, k);
        uint32_t base = kh_val(merge_from, k);

        frac->num += covg * base;
        frac->den += base;

        khint_t key = cm_put(merge_into, covg, &absent);
        if (absent) {
            kh_val(merge_into, key) = base;
        } else {
            kh_val(merge_into, key) += base;
        }
    }
}

static void process_coverage_results(covg_map *cm, fraction_t frac, char *covg_fname, char *covdist_tag, FILE *cv_file, char *cv_tag) {
    uint8_t is_nonzero = frac.num > 0 && frac.den > 0;

    double mean = -1.0;
    if (is_nonzero) {
        mean = divide(frac);
    }

    FILE *out = fopen(covg_fname, "w");
    fprintf(out, "BISCUITqc Depth Distribution - %s\ndepth\tcount\n", covdist_tag);

    khint_t k;
    uint32_t variance_numerator = 0;
    kh_foreach(cm, k) {
        uint32_t covg = kh_key(cm, k);
        uint32_t base = kh_val(cm, k);

        if (is_nonzero) {
            variance_numerator += base * (covg - mean) * (covg - mean);
        }

        fprintf(out, "%u\t%u\n", covg, base);
    }

    fflush(out);
    fclose(out);

    if (is_nonzero) {
        double sigma = sqrt((double)variance_numerator / (double)frac.den);
        fprintf(cv_file, "%s\t%lf\t%lf\t%lf\n", cv_tag, mean, sigma, sigma/mean);
    }
}

typedef struct {
    char *cv_table;     /* coefficient of variation table */
    char *all_base;     /* all base coverage */
    char *q40_base;     /* Q40 base coverage */
    char *all_cpg;      /* all cpg coverage */
    char *q40_cpg;      /* Q40 cpg coverage */
    char *all_base_top; /* all Top GC content base coverage */
    char *q40_base_top; /* Q40 Top GC content base coverage */
} output_names_t;

static inline output_names_t *init_output_names(char *prefix) {
    output_names_t *out = calloc(1, sizeof(output_names_t));

    size_t len_prefix = (prefix == NULL) ? 0 : strlen(prefix)+1;

    out->cv_table     = calloc(len_prefix + 15, sizeof(char));
    out->all_base     = calloc(len_prefix + 30, sizeof(char));
    out->q40_base     = calloc(len_prefix + 30, sizeof(char));
    out->all_cpg      = calloc(len_prefix + 30, sizeof(char));
    out->q40_cpg      = calloc(len_prefix + 30, sizeof(char));
    out->all_base_top = calloc(len_prefix + 35, sizeof(char));
    out->q40_base_top = calloc(len_prefix + 35, sizeof(char));

    if (prefix != NULL) {
        strcat(out->cv_table, prefix);
        strcat(out->all_base, prefix);
        strcat(out->q40_base, prefix);
        strcat(out->all_cpg, prefix);
        strcat(out->q40_cpg, prefix);
        strcat(out->all_base_top, prefix);
        strcat(out->q40_base_top, prefix);

        strcat(out->cv_table, "_");
        strcat(out->all_base, "_");
        strcat(out->q40_base, "_");
        strcat(out->all_cpg, "_");
        strcat(out->q40_cpg, "_");
        strcat(out->all_base_top, "_");
        strcat(out->q40_base_top, "_");
    }

    strcat(out->cv_table, "cv_table.txt");
    strcat(out->all_base, "covdist_all_base_table.txt");
    strcat(out->q40_base, "covdist_q40_base_table.txt");
    strcat(out->all_cpg, "covdist_all_cpg_table.txt");
    strcat(out->q40_cpg, "covdist_q40_cpg_table.txt");
    strcat(out->all_base_top, "covdist_all_base_topgc_table.txt");
    strcat(out->q40_base_top, "covdist_q40_base_topgc_table.txt");

    return out;
}

static inline output_names_t *destroy_output_names(output_names_t *get_wrecked) {
    free(get_wrecked->q40_base_top);
    free(get_wrecked->all_base_top);
    free(get_wrecked->q40_cpg);
    free(get_wrecked->all_cpg);
    free(get_wrecked->q40_base);
    free(get_wrecked->all_base);
    free(get_wrecked->cv_table);
    free(get_wrecked);
}

static void *coverage_write_func(void *data) {
    writer_conf_t *c = (writer_conf_t*) data;

    maps_t *maps = init_maps();

    total_coverage_t covg_fracs = init_total_coverage();
    while (1) {
        record_t rec;
        wqueue_get(record, c->q, &rec);
        if(rec.block_id == RECORD_QUEUE_END) break;

        merge(rec.maps->all_base, maps->all_base, &covg_fracs.all_base);
        merge(rec.maps->q40_base, maps->q40_base, &covg_fracs.q40_base);
        merge(rec.maps->all_cpg , maps->all_cpg , &covg_fracs.all_cpg );
        merge(rec.maps->q40_cpg , maps->q40_cpg , &covg_fracs.q40_cpg );
        if (c->has_top) {
            merge(rec.maps->all_base_top, maps->all_base_top, &covg_fracs.all_base_top);
            merge(rec.maps->q40_base_top, maps->q40_base_top, &covg_fracs.q40_base_top);
        }

        destroy_maps(rec.maps);
    }

    output_names_t *names = init_output_names(c->prefix);

    FILE *cv_table = fopen(names->cv_table, "w");
    fprintf(cv_table, "BISCUITqc Uniformity Table\ngroup\tmu\tsigma\tcv\n");

    process_coverage_results(maps->all_base, covg_fracs.all_base, names->all_base, "All Bases", cv_table, "all_base");
    process_coverage_results(maps->q40_base, covg_fracs.q40_base, names->q40_base, "Q40 Bases", cv_table, "q40_base");
    process_coverage_results(maps->all_cpg , covg_fracs.all_cpg , names->all_cpg , "All CpGs" , cv_table, "all_cpg" );
    process_coverage_results(maps->q40_cpg , covg_fracs.q40_cpg , names->q40_cpg , "Q40 CpGs" , cv_table, "q40_cpg" );
    if (c->has_top) {
        process_coverage_results(maps->all_base_top, covg_fracs.all_base_top, names->all_base_top, "All Top GC Bases", cv_table, "all_base_topgc");
        process_coverage_results(maps->q40_base_top, covg_fracs.q40_base_top, names->q40_base_top, "Q40 Top GC Bases", cv_table, "q40_base_topgc");
    }

    fflush(cv_table);
    fclose(cv_table);

    destroy_output_names(names);
    destroy_maps(maps);

    return 0;
}

static inline void increment_map(covg_map *map, khint_t bucket, int absent) {
    if (absent) {
        kh_val(map, bucket) = 1;
    } else {
        kh_val(map, bucket) += 1;
    }
}

static void format_coverage_data(maps_t *maps, uint32_t *all_covgs, uint32_t *q40_covgs, uint32_t arr_len, uint8_t *cpgs, uint8_t *tops) {
    int abs_all_base, abs_q40_base, abs_all_cpg, abs_q40_cpg, abs_all_base_top, abs_q40_base_top;
    khint_t k_all_base, k_q40_base, k_all_cpg, k_q40_cpg, k_all_base_top, k_q40_base_top;

    uint32_t i;
    for (i=0; i<arr_len; i++) {
        uint8_t is_match_all = 0;
        uint8_t is_match_q40 = 0;

        // If the current coverage is the same as the last one, we know the coverage has been seen before
        // and the correct bucket is already loaded up, so shortcircuit by automatically incrementing value
        if (i > 0) {
            is_match_all = all_covgs[i] == all_covgs[i-1];
            is_match_q40 = all_covgs[i] == all_covgs[i-1];

            if (is_match_all) {
                kh_val(maps->all_base, k_all_base) += 1;
            }
            if (is_match_q40) {
                kh_val(maps->q40_base, k_q40_base) += 1;
            }
        }

        if (!is_match_all) {
            k_all_base = cm_put(maps->all_base, all_covgs[i], &abs_all_base);
            increment_map(maps->all_base, k_all_base, abs_all_base);
        }

        if (!is_match_q40) {
            k_q40_base = cm_put(maps->q40_base, q40_covgs[i], &abs_q40_base);
            increment_map(maps->q40_base, k_q40_base, abs_q40_base);
        }

        // Always check if we're in CpG location. We'll also always pull the bucket for each hash map
        // to simplify finding buckets (rather than checking the last coverage like the all_base and
        // q40_base above)
        uint8_t is_cpg = regions_test(cpgs, i);
        if (is_cpg) {
            k_all_cpg = cm_put(maps->all_cpg, all_covgs[i], &abs_all_cpg);
            increment_map(maps->all_cpg, k_all_cpg, abs_all_cpg);

            k_q40_cpg = cm_put(maps->q40_cpg, q40_covgs[i], &abs_q40_cpg);
            increment_map(maps->q40_cpg, k_q40_cpg, abs_q40_cpg);
        }

        // Only check the top GC content if the file is provided
        // Since these regions are larger than CpGs, I may be able to keep track of the previous bin
        // and reuse that to cut down on look up time
        if (tops) {
            if (regions_test(tops, i)) {
                k_all_base_top = cm_put(maps->all_base_top, all_covgs[i], &abs_all_base_top);
                increment_map(maps->all_base_top, k_all_base_top, abs_all_base_top);

                k_q40_base_top = cm_put(maps->q40_base_top, q40_covgs[i], &abs_q40_base_top);
                increment_map(maps->q40_base_top, k_q40_base_top, abs_q40_base_top);
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

        // Extract CpGs
        uint8_t *cpgs = calloc((w.end - w.beg)/8 + 1, sizeof(uint8_t));
        regions_t *reg = get_regions(res->cpg_regions, chrm);
        if (reg) {
            int j;
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

        // Extract top GC content regions
        uint8_t *tops = calloc((w.end - w.beg)/8 + 1, sizeof(uint8_t));
        if (res->top_regions) {
            regions_t *top = get_regions(res->top_regions, chrm);
            if (top) {
                int j;
                for (j=0; j<top->n; j++) {
                    uint32_t start = top->starts[j];
                    uint32_t width = top->widths[j];
                    if (start+width >= w.beg && start < w.end) {
                        int k;
                        for (k=0; k<width; k++) {
                            if (start+k >= w.beg && start+k < w.end) {
                                regions_set(tops, start+k-w.beg);
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

            uint32_t rpos = c->pos; // 0-based reference position
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
        rec.maps = init_maps();
        format_coverage_data(rec.maps, all_covgs, q40_covgs, w.end-w.beg, cpgs, tops);

        // set record block id
        rec.block_id = w.block_id;
        // put output string to output queue
        wqueue_put2(record, res->rq, rec);

        bam_destroy1(b);
        hts_itr_destroy(iter);

        free(q40_covgs);
        free(all_covgs);
        if (tops) { free(tops); }
        free(cpgs);
    }

    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    hts_close(in);

    return 0;
}

regions_v *bed_init_regions(char *bed_fn, int n_cols) {
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
            if (strcount_char(line.s, '\t') == n_cols) {
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
    fprintf(stderr, "Usage: coverage [options] <cpgs.bed.gz> <in.bam>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -p STR    Prefix for output file names\n");
    fprintf(stderr, "    -r STR    BED file of regions of interest\n");
    fprintf(stderr, "    -s INT    Step size of windows [%d]\n", conf.step);
    fprintf(stderr, "    -@ INT    Number of threads [%d]\n", conf.n_threads);
    fprintf(stderr, "\n");

    return 1;
}

int main_coverage(int argc, char *argv[]) {
    char *prefix = 0;
    char *top_fn = 0;

    covg_conf_t conf;
    covg_conf_init(&conf);

    int c;
    if (argc < 2) { return usage(); }
    while ((c=getopt(argc, argv, ":@:p:r:s:t:")) >= 0) {
        switch (c) {
            case '@': conf.n_threads = atoi(optarg); break;
            case 'p': prefix = optarg; break;
            case 's': conf.step = atoi(optarg); break;
            case 't': top_fn = optarg; break;
            case ':': usage(); fprintf(stderr, "Option needs an argument: -%c\n", optopt); return 1;
            case '?': usage(); fprintf(stderr, "Unrecognized option: -%c\n", optopt); return 1;
            default: return usage();
        }
    }

    if (optind + 2 > argc) {
        usage();
        fprintf(stderr, "CpG BED file or BAM input is missing\n");
        return 1;
    }
    char *cpg_bed_fn = argv[optind++];
    char *infn = argv[optind++];

    regions_v *cpg_regions = bed_init_regions(cpg_bed_fn, 2);
    fprintf(stderr, "Number of chromosomes: %u\n", cpg_regions->size);

    regions_v *top_regions = top_fn ? bed_init_regions(top_fn, 3) : NULL;
    fprintf(stderr, "top_fn: %s\n", top_fn);

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
        .has_top = (top_fn) ? 1 : 0,
        .q = wqueue_init(record, conf.step),
        .prefix = prefix ? prefix : NULL,
        .targets = targets,
        .conf = &conf,
    };
    pthread_create(&writer, NULL, coverage_write_func, &writer_conf);
    for (i=0; i<conf.n_threads; ++i) {
        results[i].q = wq;
        results[i].rq = writer_conf.q;
        results[i].cpg_regions = cpg_regions;
        results[i].top_regions = top_regions;
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
        for (wbeg=0; wbeg<t->len; wbeg += conf.step, block_id++) {
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
    if (top_regions)
        destroy_regions(top_regions);

    return 0;
}
