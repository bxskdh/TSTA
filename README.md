## TSTA
pairwise and multiple sequence alignment accelerated by simd and threads
## Introduction
TSTA is a sequence alignment algorithm that combines SIMD instruction sets and thread acceleration. It uses the difference method, stripe method, and anti-diagonal method, and uses partial order alignment (POA) in multi sequence alignment.<br>

TSTA is committed to achieving faster comparison speed. TSTA only supports global alignment and allows linear and affine gap penalties. It supports SSE4.2/AVX2/AVX512 vectorization, with different choices based on hardware support.
## Installation
```
git clone https://github.com/bxskdh/TSTA.git
cd TSTA
make
```
## Usage and run
### Pairwise sequence alignment
```bash
./TSTA_psa -i ./example/psa/seq/seqa1.fa,./example/psa/seq/seqb1.fa
# maxsorce=-5
```
### no-backtracing PSA
```bash
./TSTA_psa_notrace -i ./example/psa/seq/seqa1.fa,./example/psa/seq/seqb1.fa
# maxsorce=-5
```
### Multiple sequence alignment
```bash
./TSTA_msa -i ./example/msa/seq/seq1.fa
# seq:5000, seqN:5120, poa:5000
# poa_len=5000 ,seq_len=5000 ,trace_sub:[num1]=4999 [num2]=4999 ,lastsorce=-5451
# poa_add_len:1885
# seq:5000, seqN:5120, poa:6885
# poa_len=6885 ,seq_len=5000 ,trace_sub:[num1]=6884 [num2]=4999 ,lastsorce=-3101
# poa_add_len:1714
# seq:5000, seqN:5120, poa:8599
# poa_len=8599 ,seq_len=5000 ,trace_sub:[num1]=8598 [num2]=4999 ,lastsorce=-1776
# poa_add_len:1531
# seq:5000, seqN:5120, poa:10130
# poa_len=10130 ,seq_len=5000 ,trace_sub:[num1]=10129 [num2]=4999 ,lastsorce=-870
# poa_add_len:1338
```
## Operating system
The Linux platform is supported by TSTA.It run and tested on CentOS Linux 7.
## Contact
- Peiyu Zong peiyuzong8@gmail.com (Designed and implemented core algorithms)
- Wenpeng Deng 1732889554@qq.com (Optimized some of the code) 
