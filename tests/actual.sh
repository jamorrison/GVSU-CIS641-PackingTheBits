cat all_coverage.bed | \
awk '{ cnt[$4] += $3 - $2 } END {
    for (cov in cnt) {
        sum_cov += cnt[cov]*cov
        sum_cnt += cnt[cov]
    }
    if (sum_cnt > 0 && sum_cov > 0) {
        mu = sum_cov / sum_cnt
        for (cov in cnt) {
            sum_var += cnt[cov]*((cov-mu)^2)
        }
        sigma = sqrt(sum_var / sum_cnt)
        print "all_base\t"mu"\t"sigma"\t"sigma/mu
    }
}'

cat cpg_coverage.bed | \
awk '{ cnt[$4] += 1 } END {
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
        print "all_cpg\t"mu"\t"sigma"\t"sigma/mu
    }
}'
