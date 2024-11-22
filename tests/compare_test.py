import sys

DELTA = 0.5

# Escape codes for coloring text written to terminal
GREEN='\033[0;32m'
RED='\033[0;31m'
NO_COLOR='\033[0m'

# Add green escape codes to input text
def make_green(txt):
    return f'{GREEN}{txt}{NO_COLOR}'

# Add red escape codes to input text
def make_red(txt):
    return f'{RED}{txt}{NO_COLOR}'

# Print success message
def success(txt):
    print(make_green('SUCCESS:'), txt)

# Print failure message
def fail(txt):
    print(make_red('FAIL:'), txt)

def parse_cv(fname):
    with open(fname, 'r') as fh:
        out = dict([(x.split('\t')[0], [float(n.strip()) for n in x.split('\t')[1:]]) for x in fh.readlines()[2:]])
    return out

def parse_covg(fname):
    with open(fname, 'r') as fh:
        out = dict([(int(x.split('\t')[0]), int(x.split('\t')[1])) for x in fh.readlines()[2:]])
    return out

def is_same_float(n1, n2):
    # Typical method would be to verify abs(n1-n2) < EPSILON
    # But we also need to check that two very small numbers divided by one another are the same
    # Fudge this check by going with percent difference check instead of absolute difference check
    percent_diff = abs(n1 - n2) / ((n1+n2)/2)
    if percent_diff < DELTA:
        return True

    return False

def compare_cv_tables():
    cv_orig = parse_cv(f'{sys.argv[2]}_cv_table.txt')
    cv_covg = parse_cv(f'{sys.argv[3]}_cv_table.txt')

    orig = None if 'all_base' not in cv_orig.keys() else cv_orig['all_base']
    covg = None if 'all_base' not in cv_covg.keys() else cv_covg['all_base']

    if orig is None and covg is None:
        success('All coverages are 0')
    elif orig is None and covg is not None:
        fail(f'Original was all 0 coverage, but new had non-zero coverage ({covg})')
    elif orig is not None and covg is None:
        fail(f'Original had non-zero coverage ({orig}), but new was all 0 coverage')
    else:
        if is_same_float(orig[0], covg[0]):
            success('equal means')
        else:
            fail(f'means are not equal ({orig[0]} != {covg[0]})')
        if is_same_float(orig[1], covg[1]):
            success('equal standard deviations')
        else:
            fail(f'standard deviations are not equal ({orig[1]} != {covg[1]})')
        if is_same_float(orig[2], covg[2]):
            success('equal coefficient of variations')
        else:
            fail(f'coefficient of variations are not equal ({orig[2]} != {covg[2]})')

    return None

def compare_covg_tables():
    al_orig = parse_covg(f'{sys.argv[2]}_covdist_all_base_table.txt')
    al_covg = parse_covg(f'{sys.argv[3]}_covdist_all_base_table.txt')

    for key in al_orig.keys():
        n1 = al_orig[key]
        try:
            n2 = al_covg[key]
        except KeyError:
            fail(f'coverage {key} found in original ({n1}) but not new')
            continue

        if n1 == n2:
            success(f'coverage {key} equal in both')
        else:
            fail(f'coverage {key} is not equal ({n1} != {n2})')

        al_covg.pop(key)

    if len(al_covg.keys()) > 0:
        for key in al_covg.keys():
            fail(f'coverage {key} found in new ({al_covg[key]}) but not original')

    return None

def main():
    print(f'Test: {sys.argv[1]}')
    compare_cv_tables()
    compare_covg_tables()

if __name__ == '__main__':
    main()
