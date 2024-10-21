BIN=../src/build/src

${BIN}/coverage -@ 1 -s 1000000 new.bam > tempas1.txt
valgrind --tool=memcheck --leak-check=full ${BIN}/coverage -@ 1 -s 1000000 new.bam > tempas1.txt
