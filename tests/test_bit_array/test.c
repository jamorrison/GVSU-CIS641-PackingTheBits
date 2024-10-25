#include <inttypes.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#define test(regs, i) regs[(i)>>3]&(1<<((i)&0x7))
#define set(regs, i) regs[(i)>>3] |= 1<<((i)&0x7)

void print_array(uint8_t *arr, size_t len) {
    fprintf(stdout, "[ ");

    size_t i;
    for (i=0; i<len; ++i) {
        if (i > 0) {
            fprintf(stdout, ", ");
        }
        fprintf(stdout, "%u", arr[i]);
    }

    fprintf(stdout, " ]\n");
}

void run_test(uint32_t beg, uint32_t end, uint32_t n) {
    uint32_t len = (end - beg)/8 + 1;
    fprintf(stdout, "beg: %u, end: %u, length: %u\n", beg, end, len);

    uint8_t *reg = calloc(len, sizeof(uint8_t));
    fprintf(stdout, "before: ");
    print_array(reg, (size_t)len);

    set(reg, n);
    fprintf(stdout, "after: ");
    print_array(reg, (size_t)len);

    fprintf(stdout, "n: %u, index: %u, left: %u, right: %u, logical-and: %u, test: %u\n", n, n>>3, reg[n>>3], 1<<(n&0x7), (n>>3)&(1<<(n&0x7)), test(reg, n));
}

void run_test2(uint32_t beg, uint32_t end, uint32_t lo, uint32_t hi, uint32_t bad) {
    uint32_t len = (end - beg)/8 + 1;
    fprintf(stdout, "beg: %u, end: %u, length: %u\n", beg, end, len);

    uint8_t *reg = calloc(len, sizeof(uint8_t));
    fprintf(stdout, "before: ");
    print_array(reg, (size_t)len);

    uint32_t i;
    for (i=lo; i<hi; i++) {
        set(reg, i);
        fprintf(stdout, "n = %u: ", i);
        print_array(reg, (size_t)len);

        fprintf(stdout, "n: %u, index: %u, left: %u, right: %u, logical-and: %u, test: %u\n", i, i>>3, reg[i>>3], 1<<(i&0x7), (i>>3)&(1<<(i&0x7)), test(reg, i));
    }

    fprintf(stdout, "bad test (%u): %u\n", bad, test(reg, bad));
}

void run_test3(uint32_t beg, uint32_t end, uint32_t set_start, uint32_t set_width, uint32_t *to_probe, size_t len_probe) {
    uint32_t len = (end - beg)/8 + 1;
    fprintf(stdout, "beg: %u, end: %u, length: %u\n", beg, end, len);

    uint8_t *reg = calloc(len, sizeof(uint8_t));

    uint32_t i;
    for (i=0; i<set_width; ++i) {
        set(reg, set_start + i - beg);
    }

    size_t j;
    for (j=0; j<len_probe; ++j) {
        fprintf(stdout, "\ttest value: %u, result: %u\n", to_probe[j], test(reg, to_probe[j]-beg));
    }
}

void example_01() {
    uint32_t beg = 0;
    uint32_t end = 1;

    run_test(0, 1, 0);
    run_test(0, 1, 1);
    run_test(0, 1, 2);
    run_test(0, 1, 3);
    run_test(0, 1, 4);
    run_test(0, 1, 5);
    run_test(0, 1, 6);
    run_test(0, 1, 7);

    fprintf(stdout, "\n\n");
    run_test2(0, 1, 0, 2, 7);
}

void example_02() {
    uint32_t beg = 1100000;
    uint32_t end = 1200000;
    uint32_t set = 1191893;

    uint32_t set_start = 1191893;
    uint32_t to_probe[5] = { 1191891, 1191892, 1191893, 1191894, 1191895 };

    run_test3(beg, end, set_start, 2, to_probe, 5);
    fprintf(stdout, "\n\n");
    run_test3(beg+1, end+1, set_start, 2, to_probe, 5);
}

int main() {
    //example_01();
    example_02();

    return 0;
}
