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
#include "tbx.h"

#include "khashl.h"

#include "coverage.h"

// Merge contents of one hashmap into another
// Tabulate mean coverage fraction components at the same time
static void merge(covg_map *merge_from, covg_map *merge_into, fraction_t *frac) {
    khint_t k;
    int absent;
    kh_foreach(merge_from, k) {
        uint32_t covg = kh_key(merge_from, k); // coverage level
        uint32_t base = kh_val(merge_from, k); // number of bases with coverage level

        // Tabulating weighted mean fraction components
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

// Calculate mean, standard deviation, and coefficient of variation
// Write coverage values and number of bases to output files at the same time
static void process_coverage_results(covg_map *cm, fraction_t frac, char *covg_fname, char *covdist_tag, FILE *cv_file, char *cv_tag) {
    // Only calculate statistics if we have non-zero fraction components
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

// Collect results from threads and write to various output files
static void *coverage_write_func(void *data) {
    writer_conf_t *c = (writer_conf_t*) data;

    maps_t *maps = init_maps();
    total_coverage_t covg_fracs = init_total_coverage();

    // Effectively a reduction algorithm on the hashmaps from the individual records
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
            merge(rec.maps->all_cpg_top, maps->all_cpg_top, &covg_fracs.all_cpg_top);
            merge(rec.maps->q40_cpg_top, maps->q40_cpg_top, &covg_fracs.q40_cpg_top);
        }
        if (c->has_bot) {
            merge(rec.maps->all_base_bot, maps->all_base_bot, &covg_fracs.all_base_bot);
            merge(rec.maps->q40_base_bot, maps->q40_base_bot, &covg_fracs.q40_base_bot);
            merge(rec.maps->all_cpg_bot, maps->all_cpg_bot, &covg_fracs.all_cpg_bot);
            merge(rec.maps->q40_cpg_bot, maps->q40_cpg_bot, &covg_fracs.q40_cpg_bot);
        }

        destroy_maps(rec.maps);
    }

    // Write collected results to output files
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
        process_coverage_results(maps->all_cpg_top, covg_fracs.all_cpg_top, names->all_cpg_top, "All Top GC CpGs", cv_table, "all_cpg_topgc");
        process_coverage_results(maps->q40_cpg_top, covg_fracs.q40_cpg_top, names->q40_cpg_top, "Q40 Top GC CpGs", cv_table, "q40_cpg_topgc");
    }
    if (c->has_bot) {
        process_coverage_results(maps->all_base_bot, covg_fracs.all_base_bot, names->all_base_bot, "All Bot GC Bases", cv_table, "all_base_botgc");
        process_coverage_results(maps->q40_base_bot, covg_fracs.q40_base_bot, names->q40_base_bot, "Q40 Bot GC Bases", cv_table, "q40_base_botgc");
        process_coverage_results(maps->all_cpg_bot, covg_fracs.all_cpg_bot, names->all_cpg_bot, "All Bot GC CpGs", cv_table, "all_cpg_botgc");
        process_coverage_results(maps->q40_cpg_bot, covg_fracs.q40_cpg_bot, names->q40_cpg_bot, "Q40 Bot GC CpGs", cv_table, "q40_cpg_botgc");
    }

    fflush(cv_table);
    fclose(cv_table);

    destroy_output_names(names);
    destroy_maps(maps);

    return 0;
}

// Take array of coverages and turn those coverages into a hashmap containing the coverage and number of bases
// (the number of elements with a given coverage in the array) with that coverage
static void format_coverage_data(maps_t *maps, uint32_t *all_covgs, uint32_t *q40_covgs, uint32_t arr_len, uint8_t *cpgs, uint8_t *tops, uint8_t *bots) {
    int abs_all_base, abs_q40_base, abs_all_cpg, abs_q40_cpg;
    int abs_all_base_top, abs_q40_base_top, abs_all_cpg_top, abs_q40_cpg_top;
    int abs_all_base_bot, abs_q40_base_bot, abs_all_cpg_bot, abs_q40_cpg_bot;
    khint_t k_all_base, k_q40_base, k_all_cpg, k_q40_cpg;
    khint_t k_all_base_top, k_q40_base_top, k_all_cpg_top, k_q40_cpg_top;
    khint_t k_all_base_bot, k_q40_base_bot, k_all_cpg_bot, k_q40_cpg_bot;

    uint32_t i;
    for (i=0; i<arr_len; i++) {
        uint8_t is_match_all = 0;
        uint8_t is_match_q40 = 0;

        // If the current coverage is the same as the last one, we know the coverage has been seen before
        // and the correct bucket is already in memory, so shortcircuit by automatically incrementing value
        if (i > 0) {
            is_match_all = all_covgs[i] == all_covgs[i-1];
            is_match_q40 = q40_covgs[i] == q40_covgs[i-1];

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
        uint8_t is_cpg = region_test(cpgs, i);
        if (is_cpg) {
            k_all_cpg = cm_put(maps->all_cpg, all_covgs[i], &abs_all_cpg);
            increment_map(maps->all_cpg, k_all_cpg, abs_all_cpg);

            k_q40_cpg = cm_put(maps->q40_cpg, q40_covgs[i], &abs_q40_cpg);
            increment_map(maps->q40_cpg, k_q40_cpg, abs_q40_cpg);
        }

        // Only check the top GC content if the file is provided
        //
        // Since these regions are larger than CpGs, I may be able to keep track of the previous bucket
        // and reuse that to cut down on look up time, this would require knowing when we transition
        // from not being in a region, to being in the region and then acting appropriately. I may be
        // able to set a flag when I find I'm in the region, and then turn it off when I find I'm no
        // longer in the region. If the flag is on, then check if the coverage is the same, otherwise
        // pull the bucket. This also goes for the bottom GC content regions as well.
        if (tops) {
            if (region_test(tops, i)) {
                k_all_base_top = cm_put(maps->all_base_top, all_covgs[i], &abs_all_base_top);
                increment_map(maps->all_base_top, k_all_base_top, abs_all_base_top);

                k_q40_base_top = cm_put(maps->q40_base_top, q40_covgs[i], &abs_q40_base_top);
                increment_map(maps->q40_base_top, k_q40_base_top, abs_q40_base_top);

                if (is_cpg) {
                    k_all_cpg_top = cm_put(maps->all_cpg_top, all_covgs[i], &abs_all_cpg_top);
                    increment_map(maps->all_cpg_top, k_all_cpg_top, abs_all_cpg_top);

                    k_q40_cpg_top = cm_put(maps->q40_cpg_top, q40_covgs[i], &abs_q40_cpg_top);
                    increment_map(maps->q40_cpg_top, k_q40_cpg_top, abs_q40_cpg_top);
                }
            }
        }

        // Only check the bottom GC content if the file is provided
        if (bots) {
            if (region_test(bots, i)) {
                k_all_base_bot = cm_put(maps->all_base_bot, all_covgs[i], &abs_all_base_bot);
                increment_map(maps->all_base_bot, k_all_base_bot, abs_all_base_bot);

                k_q40_base_bot = cm_put(maps->q40_base_bot, q40_covgs[i], &abs_q40_base_bot);
                increment_map(maps->q40_base_bot, k_q40_base_bot, abs_q40_base_bot);

                if (is_cpg) {
                    k_all_cpg_bot = cm_put(maps->all_cpg_bot, all_covgs[i], &abs_all_cpg_bot);
                    increment_map(maps->all_cpg_bot, k_all_cpg_bot, abs_all_cpg_bot);

                    k_q40_cpg_bot = cm_put(maps->q40_cpg_bot, q40_covgs[i], &abs_q40_cpg_bot);
                    increment_map(maps->q40_cpg_bot, k_q40_cpg_bot, abs_q40_cpg_bot);
                }
            }
        }
    }
}

static uint8_t *set_bit_array(htsFile *bed, tbx_t *tbx, int tid, uint32_t window_beg, uint32_t window_end, int n_cols) {
    // Must be free'd when done being used
    uint8_t *out = calloc((window_end - window_beg)/8 + 1, sizeof(uint8_t));
    if (bed == NULL) {
        return out;
    }

    hts_itr_t *iter = tbx_itr_queryi(tbx, tid, window_beg, window_end);
    kstring_t str = {0, 256, (char *)calloc(256, sizeof(char))};

    char *tok = NULL;
    while (tbx_itr_next(bed, tbx, iter, &str) >= 0) {
        if (strcount_char(str.s, '\t') == n_cols) {
            // Chromosome
            tok = strtok(str.s, "\t");

            // Start
            tok = strtok(NULL, "\t");
            ensure_number(tok);
            uint32_t start = (uint32_t)atoi(tok);

            // End
            tok = strtok(NULL, "\t");
            ensure_number(tok);
            uint32_t end = (uint32_t)atoi(tok);

            uint32_t width = end - start;
            uint32_t i;
            for (i=0; i<width; i++) {
                if (start+i >= window_beg && start+i < window_end) {
                    region_set(out, start + i - window_beg);
                }
            }
        }
    }

    free(str.s);
    hts_itr_destroy(iter);

    return out;
}

// Function that is run by each thread, collects coverages across genome (factoring in CIGAR string), and then
// turns the results into a hashmap for downstream processing
static void *process_func(void *data) {
    result_t *res  = (result_t*) data;

    // Open input file and check for existing BAM index
    htsFile   *in  = hts_open(res->bam_fn, "rb");
    hts_idx_t *idx = sam_index_load(in, res->bam_fn);
    if (!idx) {
        fprintf(stderr, "[%s:%d] BAM %s is not indexed?\n", __func__, __LINE__, res->bam_fn);
        fflush(stderr);
        exit(1);
    }
    bam_hdr_t *header = sam_hdr_read(in);

    record_t rec;
    memset(&rec, 0, sizeof(record_t));
    window_t w;
    while (1) {
        wqueue_get(window, res->q, &w);
        if (w.tid == -1) break;

        char *chrm = header->target_name[w.tid];


        // Tabulate coverages across window
        uint32_t *all_covgs = calloc(w.end - w.beg, sizeof(uint32_t));
        uint32_t *q40_covgs = calloc(w.end - w.beg, sizeof(uint32_t));

        // Iterate through reads that overlap window
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

            // 0-based reference position
            uint32_t rpos = c->pos;

            // Process CIGAR string to find if a base is covered or not
            int i;
            uint32_t j;
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
                        break;
                    case BAM_CINS:
                    case BAM_CSOFT_CLIP:
                        // Normally, the read position would be incremented here, but we don't actually need
                        // the read position, so leaving these separate as a historical thing rather than a
                        // practical thing
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

        // Produce coverage output
        rec.maps = init_maps();
        format_coverage_data(rec.maps, all_covgs, q40_covgs, w.end-w.beg, w.cpg, w.top, w.bot);

        // Set record block id and put output maps into output queue
        rec.block_id = w.block_id;
        wqueue_put2(record, res->rq, rec);

        // Clean up
        bam_destroy1(b);
        hts_itr_destroy(iter);

        free(q40_covgs);
        free(all_covgs);
        free(w.bot);
        free(w.top);
        free(w.cpg);
    }

    // Final clean up
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    hts_close(in);

    return 0;
}

// Print usage for tool
static int usage() {
    covg_conf_t conf;
    covg_conf_init(&conf);

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: coverage [options] <cpgs.bed.gz> <in.bam>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -p STR    Prefix for output file names\n");
    fprintf(stderr, "    -b STR    Bottom 10 percent GC content windows BED file\n");
    fprintf(stderr, "    -t STR    Top 10 percent GC content windows BED file\n");
    fprintf(stderr, "    -s INT    Step size of windows [%d]\n", conf.step);
    fprintf(stderr, "    -@ INT    Number of threads [%d]\n", conf.n_threads);
    fprintf(stderr, "    -h        Print usage\n");
    fprintf(stderr, "\n");

    return 1;
}

// Main function
int main_coverage(int argc, char *argv[]) {
    // Command line interface option values
    char *prefix = 0;
    char *top_fn = 0;
    char *bot_fn = 0;

    covg_conf_t conf;
    covg_conf_init(&conf);

    // Process command line arguments
    int c;
    if (argc < 2) { return usage(); }
    while ((c=getopt(argc, argv, ":@:b:p:s:t:h")) >= 0) {
        switch (c) {
            case '@': ensure_number(optarg); conf.n_threads = atoi(optarg); break;
            case 'b': bot_fn = optarg; break;
            case 'p': prefix = optarg; break;
            case 's': ensure_number(optarg); conf.step = atoi(optarg); break;
            case 't': top_fn = optarg; break;
            case 'h': usage(); return 0;
            case ':': usage(); fprintf(stderr, "Option needs an argument: -%c\n", optopt); return 1;
            case '?': usage(); fprintf(stderr, "Unrecognized option: -%c\n", optopt); return 1;
            default: return usage();
        }
    }

    // Missing required arguments
    if (optind + 2 > argc) {
        usage();
        fprintf(stderr, "CpG BED file or BAM input is missing\n");
        return 1;
    }

    // Set required argument values
    char *cpg_bed_fn = argv[optind++];
    char *infn = argv[optind++];

    // Validate the number of threads and step size
    if (conf.n_threads < 1) {
        fprintf(stderr, "[WARNING] Number of threads is < 1. Setting to %i\n", N_THREADS_DEFAULT);
        conf.n_threads = N_THREADS_DEFAULT;
    }
    if (conf.step < 1) {
        fprintf(stderr, "[WARNING] Window step size is < 1. Setting to %i\n", STEP_SIZE_DEFAULT);
        conf.step = STEP_SIZE_DEFAULT;
    }

    // Read input BAM to get chromosome sizes for setting windows
    htsFile *in = hts_open(infn, "rb");
    if (in == NULL) {
        fprintf(stderr, "%s unable to be opened\n", infn);
        return 1;
    }
    bam_hdr_t *header = sam_hdr_read(in);

    // Sort sequence name by alphabetic order, chr1, chr10, chr11 ...
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

    // Setup writer
    pthread_t writer;
    writer_conf_t writer_conf = {
        .has_top = (top_fn) ? 1 : 0,
        .has_bot = (bot_fn) ? 1 : 0,
        .q = wqueue_init(record, conf.step),
        .prefix = prefix ? prefix : NULL,
    };
    pthread_create(&writer, NULL, coverage_write_func, &writer_conf);

    // Setup multithreading
    wqueue_t(window) *wq = wqueue_init(window, conf.step);
    pthread_t *processors = calloc(conf.n_threads, sizeof(pthread_t));
    result_t *results = calloc(conf.n_threads, sizeof(result_t));
    for (i=0; i<conf.n_threads; ++i) {
        results[i].q = wq;
        results[i].rq = writer_conf.q;
        results[i].bam_fn = infn;
        pthread_create(&processors[i], NULL, process_func, &results[i]);
    }

    window_t w;
    memset(&w, 0, sizeof(window_t));
    uint32_t wbeg;
    int64_t block_id = 0;

    // Setup windows and add to queue for processing
    htsFile *cpg_bed = hts_open(cpg_bed_fn, "r");
    tbx_t *cpg_tbx = tbx_index_load3(cpg_bed_fn, NULL, 0);

    htsFile *top_bed = top_fn ? hts_open(top_fn, "r") : NULL;
    tbx_t *top_tbx = top_fn ? tbx_index_load3(top_fn, NULL, 0) : NULL;

    htsFile *bot_bed = bot_fn ? hts_open(bot_fn, "r") : NULL;
    tbx_t *bot_tbx = bot_fn ? tbx_index_load3(bot_fn, NULL, 0) : NULL;

    size_t j;
    for (j=0; j<targets->size; ++j) {
        t = ref_target_v(targets, j);
        for (wbeg=0; wbeg<t->len; wbeg += conf.step, block_id++) {
            w.tid = t->tid;
            w.block_id = block_id;
            w.beg = wbeg;
            w.end = wbeg + conf.step;
            if (w.end > t->len) w.end = t->len;
            w.cpg = set_bit_array(cpg_bed, cpg_tbx, t->tid, wbeg, wbeg+conf.step > t->len ? t->len : wbeg+conf.step, 2);
            w.top = set_bit_array(top_bed, top_tbx, t->tid, wbeg, wbeg+conf.step > t->len ? t->len : wbeg+conf.step, 3);
            w.bot = set_bit_array(bot_bed, bot_tbx, t->tid, wbeg, wbeg+conf.step > t->len ? t->len : wbeg+conf.step, 3);
            wqueue_put(window, wq, &w);
        }
    }

    // "Windows" for telling thread pool we've reached the end of our queue
    for (i=0; i<conf.n_threads; ++i) {
        w.tid = -1;
        wqueue_put(window, wq, &w);
    }

    // Collect our threads back into serial processing
    for (i=0; i<conf.n_threads; ++i) {
        pthread_join(processors[i], NULL);
    }

    // Tell the writer that we've reached the end of processing and can stop collecting results
    record_t rec = { .block_id = RECORD_QUEUE_END };
    wqueue_put2(record, writer_conf.q, rec);
    pthread_join(writer, NULL);

    // Clean up
    if (bot_fn) {
        tbx_destroy(bot_tbx);
        hts_close(bot_bed);
    }
    if (top_fn) {
        tbx_destroy(top_tbx);
        hts_close(top_bed);
    }
    tbx_destroy(cpg_tbx);
    hts_close(cpg_bed);
    wqueue_destroy(record, writer_conf.q);
    free_target_v(targets);
    free(results);
    free(processors);
    wqueue_destroy(window, wq);
    hts_close(in);
    bam_hdr_destroy(header);

    return 0;
}
