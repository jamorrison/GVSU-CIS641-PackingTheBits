BIN=../src/build/src
CPG=cpg.bgzip.bed.gz
TOP=windows100bp.gc_content.top10p.bgzip.bed.gz
BOT=windows100bp.gc_content.bot10p.bgzip.bed.gz
N_THREADS=1
PREFIX=example
STEP=1000000

${BIN}/coverage -b ${BOT} -t ${TOP} -p ${PREFIX} -@ "${N_THREADS}" -s ${STEP} ${CPG} example.bam
