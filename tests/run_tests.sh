BIN=../src/build/src
CPG=cpg.bgzip.bed.gz
TOP=windows100bp.gc_content.top10p.bgzip.bed.gz
BOT=windows100bp.gc_content.bot10p.bgzip.bed.gz

valgrind --tool=memcheck --leak-check=full ${BIN}/coverage -b ${BOT} -t ${TOP} -p new -@ 1 -s 1000000 ${CPG} new.bam
${BIN}/coverage -b ${BOT} -t ${TOP} -p new -@ 1 -s 1000000 ${CPG} new.bam

${BIN}/coverage -b ${BOT} -t ${TOP} -p new -@ 8 -s 1000000 ${CPG} new.bam
