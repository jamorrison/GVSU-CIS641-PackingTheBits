TIME=/varidata/research/projects/shen/tools/morrison/installed_packages/time-1.9
OUT=actual_timing.tsv
CPG=cpg.bgzip.bed.gz
#CPG=cpg.10.bed
TOP=windows100bp.gc_content.top10p.bgzip.bed.gz
BOT=windows100bp.gc_content.bot10p.bgzip.bed.gz

${TIME}/time \
    --format="genomecov_all\t%e" \
    --output=${OUT} \
    bedtools genomecov -bga -split -ibam test_files/example.bam | \
    sort -k1,1 -k2,2n > all_coverage.bed

${TIME}/time \
    --format="genomecov_q40\t%e" \
    --output=${OUT} \
    --append \
    samtools view -q 40 -b test_files/example.bam | \
    bedtools genomecov -bga -split -ibam stdin | \
    sort -k1,1 -k2,2n > q40_coverage.bed

${TIME}/time \
    --format="cpgs_all\t%e" \
    --output=${OUT} \
    --append \
    bedtools intersect \
        -sorted -wo \
        -b all_coverage.bed \
        -a ${CPG} | \
    bedtools groupby -g 1-3 -c 7 -o min > all_cpg.bed

${TIME}/time \
    --format="cpgs_q40\t%e" \
    --output=${OUT} \
    --append \
    bedtools intersect \
        -sorted -wo \
        -b q40_coverage.bed \
        -a ${CPG} | \
    bedtools groupby -g 1-3 -c 7 -o min > q40_cpg.bed

#${TIME}/time \
#    --format="cpgs_all\t%e" \
#    --output=${OUT} \
#    --append \
#    bedtools intersect \
#        -sorted -wo \
#        -b all_coverage.bed \
#        -a ${CPG} | \
#    cut -f1-3,7 > all_cpg.bed

#${TIME}/time \
#    --format="cpgs_q40\t%e" \
#    --output=${OUT} \
#    --append \
#    bedtools intersect \
#        -sorted -wo \
#        -b q40_coverage.bed \
#        -a ${CPG} | \
#    cut -f1-3,7 > q40_cpg.bed

function run_awk() {
    awk -v output=${1} -v tag=${3} '{ cnt[$4] += $3 - $2 } END {
        for (cov in cnt) {
            print int(cov)"\t"int(cnt[cov])
            sum_cov += cnt[cov]*cov
            sum_cnt += cnt[cov]
        }
        if (sum_cnt > 0 && sum_cov > 0) {
            mu = sum_cov / sum_cnt
            for (cov in cnt) {
                sum_var += cnt[cov]*((cov-mu)^2)
            }
            sigma = sqrt(sum_var / sum_cnt)
            print tag"\t"mu"\t"sigma"\t"sigma/mu >> output
        }
    }' | \
    sort -k1,1n >> ${2}
}

echo -e "BISCUITqc Uniformity Table\ngroup\tmu\tsigma\tcv" > actual_cv_table.txt

echo -e "BISCUITqc Depth Distribution - All Bases\ndepth\tcount" > actual_covdist_all_base_table.txt
${TIME}/time \
    --format="depth_all_base\t%e" \
    --output=${OUT} \
    --append \
    cat all_coverage.bed | run_awk actual_cv_table.txt actual_covdist_all_base_table.txt "all_base"

echo -e "BISCUITqc Depth Distribution - Q40 Bases\ndepth\tcount" > actual_covdist_q40_base_table.txt
${TIME}/time \
    --format="depth_q40_base\t%e" \
    --output=${OUT} \
    --append \
    cat q40_coverage.bed | run_awk actual_cv_table.txt actual_covdist_q40_base_table.txt q40_base

# Counts both bases in CpG, whereas actual BISCUITqc counts CpG as 1
echo -e "BISCUITqc Depth Distribution - All CpGs\ndepth\tcount" > actual_covdist_all_cpg_table.txt
${TIME}/time \
    --format="depth_all_cpgs\t%e" \
    --output=${OUT} \
    --append \
    cat all_cpg.bed | run_awk actual_cv_table.txt actual_covdist_all_cpg_table.txt "all_cpg"

# Counts both bases in CpG, whereas actual BISCUITqc counts CpG as 1
echo -e "BISCUITqc Depth Distribution - Q40 CpGs\ndepth\tcount" > actual_covdist_q40_cpg_table.txt
${TIME}/time \
    --format="depth_q40_cpgs\t%e" \
    --output=${OUT} \
    --append \
    cat q40_cpg.bed | run_awk actual_cv_table.txt actual_covdist_q40_cpg_table.txt "q40_cpg"

echo -e "BISCUITqc Depth Distribution - All Top GC Bases\ndepth\tcount" > actual_covdist_all_base_topgc_table.txt
${TIME}/time \
    --format="depth_all_base_topgc\t%e" \
    --output=${OUT} \
    --append \
    bedtools intersect -sorted -a all_coverage.bed -b ${TOP} | run_awk actual_cv_table.txt actual_covdist_all_base_topgc_table.txt "all_base_topgc"

echo -e "BISCUITqc Depth Distribution - Q40 Top GC Bases\ndepth\tcount" > actual_covdist_q40_base_topgc_table.txt
${TIME}/time \
    --format="depth_q40_base_topgc\t%e" \
    --output=${OUT} \
    --append \
    bedtools intersect -sorted -a q40_coverage.bed -b ${TOP} | run_awk actual_cv_table.txt actual_covdist_q40_base_topgc_table.txt "q40_base_topgc"

echo -e "BISCUITqc Depth Distribution - All Top GC CpGs\ndepth\tcount" > actual_covdist_all_cpg_topgc_table.txt
${TIME}/time \
    --format="depth_all_base_topgc\t%e" \
    --output=${OUT} \
    --append \
    bedtools intersect -sorted -a all_cpg.bed -b ${TOP} | run_awk actual_cv_table.txt actual_covdist_all_cpg_topgc_table.txt "all_cpg_topgc"

echo -e "BISCUITqc Depth Distribution - Q40 Top GC CpGs\ndepth\tcount" > actual_covdist_q40_cpg_topgc_table.txt
${TIME}/time \
    --format="depth_q40_base_topgc\t%e" \
    --output=${OUT} \
    --append \
    bedtools intersect -sorted -a q40_cpg.bed -b ${TOP} | run_awk actual_cv_table.txt actual_covdist_q40_cpg_topgc_table.txt "q40_cpg_topgc"

echo -e "BISCUITqc Depth Distribution - All Bot GC Bases\ndepth\tcount" > actual_covdist_all_base_botgc_table.txt
${TIME}/time \
    --format="depth_all_base_botgc\t%e" \
    --output=${OUT} \
    --append \
    bedtools intersect -sorted -a all_coverage.bed -b ${BOT} | run_awk actual_cv_table.txt actual_covdist_all_base_botgc_table.txt "all_base_botgc"

echo -e "BISCUITqc Depth Distribution - Q40 Bot GC Bases\ndepth\tcount" > actual_covdist_q40_base_botgc_table.txt
${TIME}/time \
    --format="depth_q40_base_botgc\t%e" \
    --output=${OUT} \
    --append \
    bedtools intersect -sorted -a q40_coverage.bed -b ${BOT} | run_awk actual_cv_table.txt actual_covdist_q40_base_botgc_table.txt "q40_base_botgc"

echo -e "BISCUITqc Depth Distribution - All Bot GC CpGs\ndepth\tcount" > actual_covdist_all_cpg_botgc_table.txt
${TIME}/time \
    --format="depth_all_base_botgc\t%e" \
    --output=${OUT} \
    --append \
    bedtools intersect -sorted -a all_cpg.bed -b ${BOT} | run_awk actual_cv_table.txt actual_covdist_all_cpg_botgc_table.txt "all_cpg_botgc"

echo -e "BISCUITqc Depth Distribution - Q40 Bot GC CpGs\ndepth\tcount" > actual_covdist_q40_cpg_botgc_table.txt
${TIME}/time \
    --format="depth_q40_base_botgc\t%e" \
    --output=${OUT} \
    --append \
    bedtools intersect -sorted -a q40_cpg.bed -b ${BOT} | run_awk actual_cv_table.txt actual_covdist_q40_cpg_botgc_table.txt "q40_cpg_botgc"

rm all_coverage.bed q40_coverage.bed all_cpg.bed q40_cpg.bed

awk '{ sum += $2 } END { print sum }' < actual_timing.tsv
