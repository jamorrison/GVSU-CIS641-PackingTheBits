function base_original() {
    test_file=${1}
    cv_table=${2}_cv_table.txt
    covg_table=${2}_covdist_all_base_table.txt
    cv_tag=${3}
    covg_tag=${4}

    echo -e "BISCUITqc Uniformity Table\ngroup\tmu\tsigma\tcv" > ${cv_table}
    echo -e "BISCUITqc Depth Distribution - ${covg_tag}\ndepth\tcount" > ${covg_table}

    bedtools genomecov -bga -split -ibam ${test_file} | \
    sort -k1,1 -k2,2n | \
    awk -v output=${cv_table} -v tag=${cv_tag} '{ cnt[$4] += $3 - $2 } END {
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
    sort -k1,1n >> ${covg_table}
}

function base_coverage() {
    test_file=${1}
    prefix=${2}
    n_threads=${3}
    step=${4}

    ../src/build/src/coverage -p ${prefix} -@ "${n_threads}" -s "${step}" cpg.bgzip.bed.gz ${test_file}
}

# Test 1 (base test)
BAM=test_files/example.bam
base_original ${BAM} a_example "all_base" "All Bases"
base_coverage ${BAM} b_example 1 100000
python compare_test.py example_test a_example b_example

# Test 2 (all_unmapped test)
BAM=test_files/all_unmapped.bam
base_original ${BAM} a_example "all_base" "All Bases"
base_coverage ${BAM} b_example 1 100000
python compare_test.py all_unmapped a_example b_example

# Test 3 (mate_unmapped test)
BAM=test_files/mate_unmapped.bam
base_original ${BAM} a_example "all_base" "All Bases"
base_coverage ${BAM} b_example 1 100000
python compare_test.py mate_unmapped a_example b_example

# Test 4 (no_flag test)
BAM=test_files/no_flag.bam
base_original ${BAM} a_example "all_base" "All Bases"
base_coverage ${BAM} b_example 1 100000
python compare_test.py no_flag a_example b_example

# Test 5 (paired test)
BAM=test_files/paired.bam
base_original ${BAM} a_example "all_base" "All Bases"
base_coverage ${BAM} b_example 1 100000
python compare_test.py paired a_example b_example

# Test 6 (proper_pair test)
BAM=test_files/paired.bam
base_original ${BAM} a_example "all_base" "All Bases"
base_coverage ${BAM} b_example 1 100000
python compare_test.py proper_pair a_example b_example

# Test 7 (dup test)
BAM=test_files/dup.bam
base_original ${BAM} a_example "all_base" "All Bases"
base_coverage ${BAM} b_example 1 100000
python compare_test.py duplicate a_example b_example

# Test 8 (qc_fail test)
BAM=test_files/qc_fail.bam
base_original ${BAM} a_example "all_base" "All Bases"
base_coverage ${BAM} b_example 1 100000
python compare_test.py qc_fail a_example b_example

# Test 9 (secondary test)
BAM=test_files/secondary.bam
base_original ${BAM} a_example "all_base" "All Bases"
base_coverage ${BAM} b_example 1 100000
python compare_test.py secondary a_example b_example

# Test 10 (supplementary test)
BAM=test_files/supplementary.bam
base_original ${BAM} a_example "all_base" "All Bases"
base_coverage ${BAM} b_example 1 100000
python compare_test.py supplementary a_example b_example
