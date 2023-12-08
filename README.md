# TSTA
pairwise and multiple sequence alignment accelerated by simd and threads
# Introduction
TSTA is a sequence alignment algorithm that combines SIMD instruction sets and thread acceleration. It uses the difference method, stripe method, and anti-diagonal method, and uses partial order alignment (POA) in multi sequence alignment.<br>

TSTA is committed to achieving faster comparison speed. TSTA only supports global alignment and allows linear and affine gap penalties. It supports SSE4.2/AVX2/AVX512 vectorization, with different choices based on hardware support.
# Installation
```
git clone https://github.com/bxskdh/TSTA.git
cd TSTA
make
```
# Usage and run
## Pairwise sequence alignment
`./TSTA_psa -M 2 -X -3 -E -2 -O -4 -T 30 -S 30 seqa.fa seqb.fa > result.txt`
## Multiple sequence alignment
`./TSTA_msa -M 2 -X -3 -E -2 -O -4 -T 30 -S 30 -i seq.fa > result.txt`
# Example
```
cd example
cd psa/msa
sh example_psa.sh/example_msa.sh
```
# Example result
## Pairwise sequence alignment
1. Result file 'r1.fa' of alignment.<br>
2. Information file 'result.txt' of alignment.
```
support:AVX512
1:10000,2:10000,length1:11520,length2:11520
B=30,T=30,maxsorce=-3
```
## Multiple sequence alignment
1. Result file 'msa.fa' of alignment.<br>
2. Consensus file 'consensus.fa' of alignment.<br>
3. Information file 'result.txt' of alignment.
```
support:AVX512
seq:5000, seqN:5760, poa:5000
poa_len=5000 ,seq_len=5000 ,trace_sub:[num1]=4999 [num2]=4999 ,lastsorce=-3378
poa_add_len:2102
seq:5000, seqN:5760, poa:7102
poa_len=7102 ,seq_len=5000 ,trace_sub:[num1]=7101 [num2]=4999 ,lastsorce=-1339
poa_add_len:1790
seq:5000, seqN:5760, poa:8892
poa_len=8892 ,seq_len=5000 ,trace_sub:[num1]=8890 [num2]=4999 ,lastsorce=-98
poa_add_len:1520
seq:5000, seqN:5760, poa:10412
poa_len=10412 ,seq_len=5000 ,trace_sub:[num1]=10410 [num2]=4999 ,lastsorce=766
poa_add_len:1341
```
# Operating system
The Linux platform is supported by TSTA.It run and tested on CentOS Linux 7.
# Contact
Peiyu Zong peiyuzong8@gmail.com
