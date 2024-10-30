import subprocess
import shutil
import os

# Files / executable
COVG = '../../src/build/src/coverage'
CPG  = '../cpg.bgzip.bed.gz'
BAM  = '../example.bam'

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

# Create directory and move into it
def make_dir_and_move(dir_name):
    os.mkdir(dir_name)
    os.chdir(dir_name)

# Move to parent directory and delete input directory
def move_up_and_delete(dir_name):
    os.chdir('../')
    shutil.rmtree(dir_name)

def run_test_01():
    """Check for optional command line parameter for number of threads."""
    test_dir = 'test_01'
    make_dir_and_move(test_dir)

    cmd1 = f'{COVG} -p test_01a {CPG} {BAM}'         # default
    cmd2 = f'{COVG} -p test_01b -@ 1 {CPG} {BAM}'    # include argument
    cmd3 = f'{COVG} -p test_01c -@ test {CPG} {BAM}' # bad input

    res1 = subprocess.run(cmd1, capture_output=True, text=True, shell=True)
    res2 = subprocess.run(cmd2, capture_output=True, text=True, shell=True)
    res3 = subprocess.run(cmd3, capture_output=True, text=True, shell=True)

    if res1.returncode == 0 and res2.returncode == 0 and res3.returncode == 1:
        success('Test 01')
    else:
        fail('Test 01')

    move_up_and_delete(test_dir)

def main():
    run_test_01()

if __name__ == '__main__':
    main()
