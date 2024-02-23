CC = gcc
OPTIMIZE_FLAGS = -O3 -g
SIMD_FLAGS = -msse4.2 #optional:-march=native
THREAD_FLAGS = -lpthread

THREAD_SITE = pthreadpool
THREAD_SRC = ${THREAD_SITE}/pthreadpool.h
OTHER_SRC = msa/poa.h

all : TSTA_psa TSTA_psa_notrace TSTA_msa clean

#Pairwise sequence alignment

TSTA_psa : ${THREAD_SRC} psa.o pthreadpool.o seqio.o
	${CC} ${OPTIMIZE_FLAGS} -o $@ $^ -lz -lm ${THREAD_FLAGS}

psa.o : psa/psa.c
	${CC} ${OPTIMIZE_FLAGS} -c psa/psa.c -o $@ ${SIMD_FLAGS}

#psa-notrace

TSTA_psa_notrace : ${THREAD_SRC} psa_notrace.o pthreadpool.o
	${CC} ${OPTIMIZE_FLAGS} -o $@ psa_notrace.o pthreadpool.o ${THREAD_FLAGS}

psa_notrace.o: psa/psa_notrace.c
	${CC} ${OPTIMIZE_FLAGS} -c psa/psa_notrace.c -o $@ ${SIMD_FLAGS}

#Multiple sequence alignment

TSTA_msa : ${THREAD_SRC} ${OTHER_SRC} msa.o c-t-simd.o topo.o result.o pthreadpool.o seqio.o
	${CC} ${OPTIMIZE_FLAGS} -o $@ msa.o c-t-simd.o topo.o result.o pthreadpool.o seqio.o -lz -lm ${THREAD_FLAGS}

msa.o : msa/msa.c ${OTHER_SRC}
	${CC} ${OPTIMIZE_FLAGS} -c msa/msa.c -o $@ ${SIMD_FLAGS}

c-t-simd.o : msa/c-t-simd.c ${OTHER_SRC}
	${CC} ${OPTIMIZE_FLAGS} -c msa/c-t-simd.c -o $@ ${SIMD_FLAGS}

topo.o : msa/topo.c ${OTHER_SRC}
	${CC} ${OPTIMIZE_FLAGS} -c msa/topo.c -o $@

result.o : msa/result.c ${OTHER_SRC}
	${CC} ${OPTIMIZE_FLAGS} -c msa/result.c -o $@

#pthreadpool

pthreadpool.o : ${THREAD_SITE}/pthreadpool.c
	${CC} ${OPTIMIZE_FLAGS} -c ${THREAD_SITE}/pthreadpool.c -o $@ ${THREAD_FLAGS}

seqio.o: seqio.c
	${CC} ${OPTIMIZE_FLAGS} -c seqio.c -lz -lm -o $@

#clean
clean:
	rm -f *.o

