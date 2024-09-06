#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "poa.h"
#include "../pthreadpool/pthreadpool.h"
#include <unistd.h>
#include <semaphore.h>
#include <pthread.h>
#include <limits.h>
#include <malloc.h>
#include <stdint.h>

#define nconvert(n) (maxtag > 0 ? (!!(n / maxtag)) * maxtag + (!(n / maxtag)) * (n % maxtag) : 0)
#define NUM2(num2) ((num2) / L) * L + ((((num2) % L) % W) * B + (((num2) % L) / W))

int s_len;//temp
char z = 0;
extern char M, X, E, O;
extern int bS, L, B, W;

int maxtag, fmaxtag, lmaxtag, length1, length2;
volatile int lock;

pthread_mutex_t mutex;

typedef struct pack
{
	int i;
	int j;
	int l;
	int num;
	char* seq;
	topo* p;
}pack;

char* readseq(char* A, int poal)
{
	char* seq1, * seq2;
	int a = strlen(A);
	length1 = a;
	length2 = poal;
	if (a % L != 0)
		length1 = a + (L - a % L);
	seq1 = (char*)malloc(((size_t)length1 + 1) * sizeof(char));
	seq2 = (char*)mm_malloc((size_t)(length1+1) * sizeof(char));
	memset(seq1, 'N', (size_t)length1 + 1);
	memcpy(seq1, A, a);
	seq1[length1] = '\0';
	for (int i = 0; i < length1; i++)
		seq2[i] = seq1[i / L * L + ((i % L) % B) * W + ((i % L) / B)];
	free(seq1);
	seq2[length1] = '\0';
	return seq2;
}

void free_node(poa* n)
{
	if (n->pre)
	{
		free(n->pre);
		n->pre = NULL;
	}
	if (n->next)
	{
		free(n->next);
		n->next = NULL;
	}
	free(n);
	n = NULL;
}

poa* poa_build_init(topo* p, char a[], int sum)
{
	int len_a = strlen(a);
	p->len = len_a;
	p->sort = (poa**)malloc(len_a * sizeof(poa*));
	p->unsort = (poa**)malloc(len_a * sizeof(poa*));
	p->last_node_num = 1;
	poa* head, * node, * end, * init;

	init = (poa*)malloc(sizeof(poa));
	init->frist_col_sorce = 0;
	init->sorce = (char*)mm_malloc(L * sizeof(char));
	init->esorce = (char*)mm_malloc(L * sizeof(char));
	memset(init->sorce, E, L);//u
	memset(init->esorce, E + E + O, L);//e
	init->simple_sorce = (int*)malloc(sizeof(int));
	init->simple_sorce[0] = 0;
	init->sub = init->node_logo = -1;
	init->in = -1;
	init->out = -1;
	init->base = 'N';
	init->node_sorce = init->node_base_len = 0;
	init->pre = NULL; init->next = NULL; init->source = NULL; init->edge_weight = NULL; init->f0 = NULL;

	head = (poa*)malloc(sizeof(poa));
	head->pre = (poa**)malloc(sizeof(poa*));
	head->next = (poa**)malloc(sizeof(poa*));
	head->next[0] = (poa*)malloc(sizeof(poa));
	head->source = head->esource = head->fsource = NULL;
	head->sorce = (char*)mm_malloc(L * sizeof(char));
	head->esorce = (char*)mm_malloc(L * sizeof(char));
	head->simple_sorce = NULL;
	head->pre[0] = init;//
	head->base = a[0];
	head->passing_seq = (char*)malloc(sum * sizeof(char));
	memset(head->passing_seq, 0, sum * sizeof(char));
	head->passing_seq[0] = 1;
	head->in = 0;
	head->out = 1;
	head->sub = 0;
	head->node_logo = 0;
	head->node_sorce = 0;
	head->node_base_len = 1;
	head->edge_weight = (int*)malloc(sizeof(int));
	head->edge_weight[0] = 0;
	head->f0 = (char*)malloc(sizeof(char));
	head->mismatch_num = 0;
	p->unsort[0] = p->sort[0] = head;

	end = head;
	for (int i = 1; i < len_a; i++)
	{
		node = end->next[0];
		if (node)
		{
			node->pre = (poa**)malloc(sizeof(poa*));
			node->pre[0] = end;
			node->next = (poa**)malloc(sizeof(poa*));
			node->next[0] = (poa*)malloc(sizeof(poa));
			node->source = node->esource = node->fsource = NULL;
			node->sorce = (char*)mm_malloc(L * sizeof(char));
			node->esorce = (char*)mm_malloc(L * sizeof(char));
			node->simple_sorce = NULL;
			node->base = a[i];
			node->in = 1;
			node->out = 1;
			node->sub = i;
			node->passing_seq = (char*)malloc(sum * sizeof(char));
			memset(node->passing_seq, 0, sum * sizeof(char));
			node->passing_seq[0] = 1;
			node->node_logo = node->node_sorce = 0;
			node->edge_weight = (int*)malloc(sizeof(int));
			node->edge_weight[0] = 1;
			node->f0 = (char*)malloc(sizeof(char));
			node->mismatch_num = 0;
			p->unsort[i] = p->sort[i] = node;

			end = node;
		}
		else
		{
			printf("init wrong!");
			exit(0);
		}
	}

	free(end->next);
	end->next = NULL;
	end->out = 0;

	return head;
}

void block_line_alignment(int block_i, int block_j, int block_l, poa* row, char** f_temp, char** r_s, char* h_g, char* v0, char* seq, int nv, int pc2, int* pd, int* te, char** VC2, char** VC1, char* vc_1, char* vc_2)
{
	int m1, m2, m3;
	m1 = m2 = m3 = 0;
	char logo = -6;
	char Logo1 = 60;
	char Logo = 100;
	short reduce = 0;
	int pre_num = row->in;
	if (pre_num == 0)
	{
		if (block_i == 0)
		{
			row->pre[0]->sorce[0] = O + E;
			row->pre[0]->esorce[0] = 2 * (O + E);
		}
		else
		{
			row->pre[0]->sorce[0] = E;
			row->pre[0]->esorce[0] = E + O + E;
		}
		pre_num = 1;
		row->frist_col_sorce = row->simple_sorce[0] = E + O;
	}
		
	for (int i = 0; i < pre_num; i++)
		pd[i] = (row->pre[i]->node_logo / 3) * pc2;
	int pc1 = (row->node_logo / 3) * pc2;

	if (block_i <= lmaxtag && block_l == block_j - 1 && row->in != 0)
	{
		row->frist_col_sorce = row->pre[0]->frist_col_sorce + E;
		for (int i = 1; i < pre_num; i++)
			if (row->frist_col_sorce < row->pre[i]->frist_col_sorce + E)
				row->frist_col_sorce = row->pre[i]->frist_col_sorce + E;
		row->simple_sorce[0] = row->frist_col_sorce;
		for (int i = 0; i < pre_num; i++)
		{
			te[i] = row->frist_col_sorce - row->pre[i]->frist_col_sorce;
			if (te[i] > Logo)
			{
				v0[i] = Logo;
				if(te[i] - Logo > 127)
				{
					vc_2[i] = VC2[i][0] = (te[i] - Logo - 127) > 127 ? 127 : (te[i] - Logo - 127);
					vc_1[i] = VC1[i][0] = 127;
				}
				else
				{
					vc_2[i] = VC2[i][0] = 0;
					vc_1[i] = VC1[i][0] = te[i] - Logo;
				}
			}
			else
			{
				v0[i] = te[i];
				vc_2[i] = VC2[i][0] = 0;
				vc_1[i] = VC1[i][0] = 0;
			}
		}
	}
	else
	{
		if (row->pre[0]->sub == -1)
		{
			v0[0] = row->simple_sorce[nv] - (nv * L * E + (nv > 0 ? O : 0));
			vc_2[0] = VC2[0][0] = 0;
			vc_1[0] = VC1[0][0] = 0;
		}
		else
		{
			for (int i = 0; i < pre_num; i++)
			{
				te[i] = row->simple_sorce[nv] - row->pre[i]->simple_sorce[nv];
				if (te[i] > Logo)
				{
                                       	v0[i] = Logo;
	                                if(te[i] - Logo > 127)
        	                        {
                	                        vc_2[i] = VC2[i][0] = (te[i] - Logo - 127) > 127 ? 127 : (te[i] - Logo - 127);
                        	                vc_1[i] = VC1[i][0] = 127;
                        	        }
                                	else
                                	{
                                        	vc_2[i] = VC2[i][0] = 0;
                                        	vc_1[i] = VC1[i][0] = te[i] - Logo;
                                	}
				}
				else
				{
					v0[i] = te[i];
					vc_2[i] = VC2[i][0] = 0;
                               		vc_1[i] = VC1[i][0] = 0;
				}
			}
		}	
	}

	if (block_i <= lmaxtag && block_l == block_j - 1 && block_i < length2 / L)
	{
		if(row->in == 0)
			row->f0[0] = v0[0] + E + O;
		else
			for (int i = 0; i < row->in; i++)
				row->f0[i] = v0[i] + E + O;
	}
	__mxxxi zero, top, Smin, h, max, emax, eumax, source, source_num, esource, esource_num, fsource, s0, s1, s2, s3, s4, s5, s6, temp, temp1, temp2, mat, mis, egap, ogap, base, N, z;
	__mask mask, mask1, mask2, mask3, mask4, mask5, SN ,SM, SX;
	__mxxxi y[pre_num], vc2[pre_num], vc1[pre_num], vc0[pre_num], diff[pre_num], t[pre_num], e[pre_num], eu[pre_num], ev[pre_num], f[pre_num], fv[pre_num], q[pre_num], v[pre_num], sum[pre_num];
	z = mm_set1_epi8(Logo1);
	zero = mm_setzero();
	top = mm_set1_epi8(127);
	Smin = mm_set1_epi8(MIN);	
	for (int i = 0; i < pre_num; i++)
	{
		sum[i] = mm_setzero();
		for (int j = 0; j < W; j++)
		{
			s1 = mm_load((__mxxxi*)row->pre[i]->sorce + pd[i] + j);
			sum[i] = mm_add_epi8(sum[i], s1);
		}
		mm_store((__mxxxi*)r_s[i], sum[i]);
	}
        
	if (pre_num != 1)
        {
                for (int i = 0; i < pre_num; i++)
                        f_temp[i][0] = v0[i];
                for (int j = 1; j < B; j++)
                {
			for (int i = 0; i < pre_num; i++)
			{
				te[i] = te[i] - r_s[i][j - 1] + W * E;
			}
			m1 = te[0];
			for (int s = 1; s < pre_num; s++)
			{
				if (te[s] < m1)
					m1 = te[s];
			}
			m2 = logo - m1;
			for (int i = 0; i < pre_num; i++)
			{
				if (te[i] + m2 > Logo)
				{
					f_temp[i][j] = Logo;
					if(te[i] + m2 - Logo > 127)
                                        {
                                                VC2[i][j] = (te[i] + m2 - Logo - 127) > 127 ? 127 : (te[i] + m2 - Logo - 127);
                                                VC1[i][j] = 127;
                                        }
                                        else
                                        {
                                                VC2[i][j] = 0;
                                                VC1[i][j] = te[i] + m2 - Logo;
                                        }
				}
				else
				{
					f_temp[i][j] = te[i] + m2;
					VC2[i][j] = 0;
					VC1[i][j] = 0;
				}
			}
                }
                for (int i = 0; i < pre_num; i++)
                        v[i] = mm_load((__mxxxi*)f_temp[i]);
        }
	else
	{
		for (int j = 0; j < B; j++)
		{
			VC2[0][j] = 0;
			VC1[0][j] = 0;
		}
		vc_1[0] = vc_2[0] = 0;
		v[0] = mm_set1_epi8(E);
		v[0] = mm_insert_epi8(v[0], v0[0], 0);
	}
	mat = mm_set1_epi8(M);
	mis = mm_set1_epi8(X);
	egap = mm_set1_epi8(E);
	ogap = mm_set1_epi8(O + E);
	base = mm_set1_epi8(row->base);
	N = mm_set1_epi8('N');
	for (int j = 0; j < pre_num; j++)
	{
		vc2[j] = mm_load((__mxxxi*)VC2[j]);
		vc1[j] = mm_load((__mxxxi*)VC1[j]);///
		f[j] = Smin;
		f[j] = mm_insert_epi8(f[j], row->f0[j], 0);
	}
	for (int i = 0; i < W; i++)
	{
		h = mm_load((__mxxxi*)seq + pc2 + i);
		mask = mm_cmpeq_epi8(h, base);
		h = mm_blendv_epi8(mis, mat, mask);
		mm_store((__mxxxi*)h_g + i, h);
		s1 = Smin;
		for (int j = 0; j < pre_num; j++)
		{
			t[j] = mm_load((__mxxxi*)row->pre[j]->sorce + pd[j] + i);
			e[j] = mm_load((__mxxxi*)row->pre[j]->esorce + pd[j] + i);
			temp = mm_max_epi8(f[j], h);
			temp = mm_max_epi8(e[j], temp);
			temp = mm_subs_epi8(temp, v[j]);
			mask4 = mm_cmpgt_epi8(v[j], z);
			temp = mm_blendv_epi8(temp, ogap, mask4);
			s1 = mm_max_epi8(s1, temp);
		}
		for (int j = 0; j < pre_num; j++)
		{
			temp = mm_sub_epi8(t[j], egap);
			temp = mm_subs_epi8(f[j], temp);
			temp1 = mm_adds_epi8(s1, ogap);
			temp1 = mm_subs_epi8(temp1, t[j]);
			temp1 = mm_adds_epi8(v[j], temp1);
			f[j] = mm_max_epi8(temp, temp1);

			temp1 = mm_subs_epi8(s1, t[j]);
			vc0[j] = mm_adds_epi8(v[j], temp1);

			mask4 = mm_cmpgt_epi8(temp1, zero);
			temp1 = mm_blendv_epi8(zero, temp1, mask4);
			temp2 = mm_subs_epi8(top, v[j]);
			y[j] = mm_subs_epu8(temp1, temp2);

			v[j] = mm_adds_epi8(vc0[j], vc1[j]);
			
			mask5 = mm_cmpeq_epi8(vc1[j], zero);
			temp2 = mm_subs_epu8(top, vc0[j]);
			diff[j] = mm_blendv_epi8(temp2, zero, mask5);
			
			temp2 = vc1[j];
			vc1[j] = mm_subs_epu8(vc1[j], diff[j]);
			vc1[j] = mm_adds_epi8(vc1[j], vc2[j]); 
			temp2 = mm_subs_epu8(vc1[j], temp2);
			vc2[j] = mm_subs_epu8(vc2[j], diff[j]);
			vc2[j] = mm_adds_epi8(vc2[j], y[j]);
			vc2[j] = mm_subs_epu8(vc2[j], temp2);
		}
	}	

	for (int j = 0; j < pre_num; j++)
	{
		mm_store((__mxxxi*)f_temp[j], f[j]);
		te[j] = f_temp[j][0];
		for (int x = 1; x < B - 1; x++)
		{
			te[j] = te[j] - r_s[j][x] + W * E;
			if (te[j] > f_temp[j][x] && te[j] > 125)
			{
				f_temp[j][x] = 125;
			}
			else if (te[j] > f_temp[j][x] && te[j] <= 125)
			{
				f_temp[j][x] = te[j];
			}
			else if (te[j] <= f_temp[j][x] && f_temp[j][x] > 125)
			{
				te[j] = f_temp[j][x];
				f_temp[j][x] = 125;
			}
			else
			{
				te[j] = f_temp[j][x];
			}
		}
		f[j] = mm_load((__mxxxi*)f_temp[j]);
		temp1 = mm_subs_epi8(f[j], egap);
		f[j] = mm_slli(f[j]);
		f[j] = mm_insert_epi8(f[j], row->f0[j], 0);

		vc0[j] = mm_max_epi8(temp1, v[j]);
		vc0[j] = mm_slli(vc0[j]);
		vc0[j] = mm_insert_epi8(vc0[j], v0[j], 0);
		
                vc1[j] = mm_slli(vc1[j]);
                vc1[j] = mm_insert_epi8(vc1[j], vc_1[j], 0);
		v[j] = mm_adds_epi8(vc0[j], vc1[j]);

		vc2[j] = mm_slli(vc2[j]);
                vc2[j] = mm_insert_epi8(vc2[j], vc_2[j], 0);
	}
	
	sum[0] = zero;
	s0 = mm_set1_epi8(42);
	s2 = mm_add_epi8(s0, s0);
	s3 = mm_add_epi8(s2, s0);
	s4 = mm_set1_epi8(1);
	s5 = mm_add_epi8(s0, s4);
	s6 = mm_add_epi8(s4, s4);

	for (int i = 0; i < W; i++)
	{
		h = mm_load((__mxxxi*)seq + pc2 + i);
		SN = mm_cmpeq_epi8(h, N);
		h = mm_load((__mxxxi*)h_g + i);
		SM = mm_cmpeq_epi8(mat, h);
		SX = mm_cmpeq_epi8(mis, h);
		max = eumax = Smin;
		for (int j = 0; j < pre_num; j++)
		{
			t[j] = mm_load((__mxxxi*)row->pre[j]->sorce + pd[j] + i);
			e[j] = mm_load((__mxxxi*)row->pre[j]->esorce + pd[j] + i);
			fv[j] = mm_subs_epi8(f[j], v[j]);
			eu[j] = mm_subs_epi8(e[j], v[j]);
			q[j] = mm_subs_epi8(h, v[j]);
			temp = mm_max_epi8(fv[j], eu[j]);
			temp = mm_max_epi8(temp, q[j]);
			mask4 = mm_cmpgt_epi8(v[j], z);
                        temp = mm_blendv_epi8(temp, ogap, mask4);
			max = mm_max_epi8(max, temp);
			ev[j] = mm_subs_epi8(e[j], t[j]);
			eumax = mm_max_epi8(eumax, eu[j]);
		}
		max = mm_blendv_epi8(max, zero, SN);
		sum[0] = mm_add_epi8(sum[0], max);

		//source
		source = s3;
		source_num = zero;
		for (int j = pre_num - 1; j >= 0; j--)
		{
			mask = mm_cmpeq_epi8(max, eu[j]);
			source = mm_blendv_epi8(source, zero, mask);
			source_num = mm_blendv_epi8(source_num, mm_set1_epi8(j), mask);
		}
		for (int j = pre_num - 1; j >= 0; j--)
		{
			mask = mm_and_epi8(mm_cmpeq_epi8(max, q[j]), SX);
			source = mm_blendv_epi8(source, s2, mask);
			source_num = mm_blendv_epi8(source_num, mm_set1_epi8(j), mask);
		}
		for (int j = pre_num - 1; j >= 0; j--)
		{
			mask = mm_and_epi8(mm_cmpeq_epi8(max, q[j]), SM);
			source = mm_blendv_epi8(source, s0, mask);
			source_num = mm_blendv_epi8(source_num, mm_set1_epi8(j), mask);
		}
		source = mm_add_epi8(source, source_num);//now:0/42/84/126
		mm_store((__mxxxi*)row->source + pc2 + i, source);
		mm_store((__mxxxi*)row->sorce + pc1 + i, max);
		
		//esource+fsource
		esource = fsource = s4;
		esource_num = zero;
		temp = mm_adds_epi8(max, ogap);
		emax = Smin;
		for (int j = pre_num - 1; j >= 0; j--)
		{
			f[j] = mm_adds_epi8(f[j], egap);
			s1 = mm_adds_epi8(temp, v[j]);
			mask1 = mm_cmpeq_epi8(f[j], s1);
			f[j] = mm_max_epi8(f[j], s1);
			f[j] = mm_subs_epi8(f[j], t[j]);
			mask = mm_cmpeq_epi8(fv[j], ogap);
			fsource = mm_blendv_epi8(fsource, s6, mask);

			e[j] = mm_adds_epi8(e[j], egap);
			e[j] = mm_subs_epi8(e[j], v[j]);
			mask2 = mm_cmpeq_epi8(temp, e[j]);
			temp1 = mm_max_epi8(temp, e[j]);
			emax = mm_max_epi8(emax, temp1);

			mask3 = mm_cmpeq_epi8(eu[j], eumax);
			esource_num = mm_blendv_epi8(esource_num, mm_set1_epi8(j), mask3);
			mask = mm_cmpeq_epi8(ev[j], ogap);
			temp1 = mm_blendv_epi8(s4, s5, mask);
			temp1 = mm_add_epi8(temp1, esource_num);
			esource = mm_blendv_epi8(esource, temp1, mask3);
			temp1 = mm_sub_epi8(zero, esource);
			mask = mm_and_epi8(mask3, mask2);
			esource = mm_blendv_epi8(esource, temp1, mask);

                        temp1 = mm_subs_epi8(max, t[j]);
                        vc0[j] = mm_adds_epi8(v[j], temp1);

			mask4 = mm_cmpgt_epi8(temp1, zero);
			temp1 = mm_blendv_epi8(zero, temp1, mask4);
                        temp2 = mm_subs_epi8(top, v[j]);
                        y[j] = mm_subs_epu8(temp1, temp2);

                        v[j] = mm_adds_epi8(vc0[j], vc1[j]);

                        mask5 = mm_cmpeq_epi8(vc1[j], zero);
                        temp2 = mm_subs_epu8(top, vc0[j]);
                        diff[j] = mm_blendv_epi8(temp2, zero, mask5);

                        temp2 = vc1[j];
                        vc1[j] = mm_subs_epu8(vc1[j], diff[j]);
                        vc1[j] = mm_adds_epi8(vc1[j], vc2[j]); 
                        temp2 = mm_subs_epu8(vc1[j], temp2);
                        vc2[j] = mm_subs_epu8(vc2[j], diff[j]);
                        vc2[j] = mm_adds_epi8(vc2[j], y[j]);
                        vc2[j] = mm_subs_epu8(vc2[j], temp2);
			/*temp1 = mm_subs_epi8(max, t[j]);
                        v[j] = mm_adds_epi8(v[j], temp1);*/
		}
		temp1 = mm_sub_epi8(zero, fsource);
		fsource = mm_blendv_epi8(fsource, temp1, mask1);
		mm_store((__mxxxi*)row->fsource + pc2 + i, fsource);//-
		mm_store((__mxxxi*)row->esource + pc2 + i, esource);//-
		mm_store((__mxxxi*)row->esorce + pc1 + i, emax);//-
	}
	for (int j = 0; j < pre_num; j++)
		row->f0[j] = mm_extract_epi8(f[j]);
	s1 = mm_reduceto16_epi8(sum[0]);
	s1 = mm_hadd_epi8(s1);
	reduce = mm_reduce_epi8(s1);
	row->simple_sorce[nv + 1] = row->simple_sorce[nv] + reduce;

        if (row->out == 0 && block_i >= maxtag && block_l == 0)
                row->lastsorce = row->simple_sorce[nv + 1];

	//cross-block
	int kk = (row->sub / L + 1) * L;
	for (int i = 0; i < row->out; i++)
	{
		if (row->next[i]->sub >= kk && row->node_logo != 3)
		{
			char* t_sorce = (char*)mm_malloc(length1 * sizeof(char));//
                        memcpy(t_sorce, row->sorce, L * sizeof(char));
                        mm_free(row->sorce);
                        row->sorce = t_sorce;

			char* t_esorce = (char*)mm_malloc(length1 * sizeof(char));//
			memcpy(t_esorce, row->esorce, L * sizeof(char));
			mm_free(row->esorce);
			row->esorce = t_esorce;

                        row->node_logo = 3;
		}
	}
}

//void block_alignment(int block_i, int block_j, int block_l, topo* p, char* seq)
void block_alignment(void* pa)
{
	pack* pb = (pack*)pa;
	int block_i = pb->i;
	int block_j = pb->j;
	int block_l = pb->l;
	int num = pb->num;
	char* seq = pb->seq;
	topo* p = pb->p;

	int a1, a2;
	int nv = nconvert(block_i) - block_l;
	int pc2 = nv * L / B;

	char** f_temp = (char**)malloc(num * sizeof(char*));///
        for (int i = 0; i < num; i++)
                f_temp[i] = (char*)mm_malloc(B * sizeof(char));

        char** VC2 = (char**)malloc(num * sizeof(char*));///
        for (int i = 0; i < num; i++)
                VC2[i] = (char*)mm_malloc(B * sizeof(char));
	
	char** VC1 = (char**)malloc(num * sizeof(char*));///
        for (int i = 0; i < num; i++)
                VC1[i] = (char*)mm_malloc(B * sizeof(char));

	char** r_s = (char**)malloc(num * sizeof(char*));///g
	for (int i = 0; i < num; i++)
		r_s[i] = (char*)mm_malloc(B * sizeof(char));
	char* h_g = (char*)mm_malloc(L * sizeof(char));
	char* v0 = (char*)malloc(num * sizeof(char));
	char* vc_1 = (char*)malloc(num * sizeof(char));
	char* vc_2 = (char*)malloc(num * sizeof(char));
	int* pd = (int*)malloc(num * sizeof(int));
	int* te = (int*)malloc(num * sizeof(int));
	a1 = (((block_i - maxtag) > 0) * (block_i - maxtag) + block_l) * L;
	for (int i = 0; i < L; i++)
	{
		a2 = a1 + i;
		if (a2 >= p->len)
			break;
		block_line_alignment(block_i, block_j, block_l, p->sort[a2], f_temp, r_s, h_g, v0, seq, nv, pc2, pd, te, VC2, VC1, vc_1, vc_2);
	}
	for (int i = 0; i < num; i++)
	{
		mm_free(VC2[i]);
		mm_free(VC1[i]);
		mm_free(f_temp[i]);
		mm_free(r_s[i]);
	}
	free(VC2);free(VC1);free(f_temp);free(r_s);
	free(v0); mm_free(h_g); free(pd); free(te); free(vc_1); free(vc_2);
        pthread_mutex_lock(&mutex);
        lock++;
        pthread_mutex_unlock(&mutex);
}

topo* node_fuse(topo* n, char b[], int num, int sum, int last)
{
	int s3 = 0;
	int l = 0;
	poa* init;
	init = (poa*)malloc(sizeof(poa));
	init->pre = (poa**)malloc(sizeof(poa*));
	init->frist_col_sorce = 0;
	init->source = NULL;
	init->next = NULL;
	init->sorce = (char*)mm_malloc(L * sizeof(char));
	init->esorce = (char*)mm_malloc(L * sizeof(char));
	memset(init->sorce, E, L);
	memset(init->esorce, E + E + O, L);
	init->simple_sorce = (int*)malloc(sizeof(int));
	init->simple_sorce[0] = 0;
	init->sub = init->node_logo = -1;
	init->pre[0] = NULL;
	init->in = -1;
	init->node_sorce = init->node_base_len = 0;
	init->edge_weight = NULL;
	init->f0 = NULL;

	int len_b = strlen(b);
	n->unsort = (poa**)realloc(n->unsort, (len_b + n->len) * sizeof(poa*));

	poa** seq = (poa**)malloc(len_b * sizeof(poa*));
	seq[0] = (poa*)malloc(sizeof(poa));
	seq[0]->pre = (poa**)malloc(sizeof(poa*));
	seq[0]->pre[0] = init;
	seq[0]->next = (poa**)malloc(sizeof(poa*));
	seq[0]->base = b[0];
	seq[0]->source = seq[0]->esource = seq[0]->fsource = NULL;
	seq[0]->sorce = seq[0]->esorce = NULL;
	seq[0]->simple_sorce = NULL;
	seq[0]->node_logo = 0;
	seq[0]->in = 0;
	seq[0]->out = 1;
	seq[0]->sub = -1;
	seq[0]->node_sorce = 0;
	seq[0]->node_base_len = 1;
	seq[0]->edge_weight = NULL;
	seq[0]->passing_seq = NULL;
	seq[0]->mismatch_num = 0;

	for (int i = 1; i < len_b; i++)
	{
		seq[i] = (poa*)malloc(sizeof(poa));
		seq[i]->pre = (poa**)malloc(sizeof(poa*));
		seq[i]->pre[0] = seq[i - 1];
		seq[i - 1]->next[0] = seq[i];
		seq[i]->next = (poa**)malloc(sizeof(poa*));
		seq[i]->base = b[i];
		seq[i]->source = seq[i]->esource = seq[i]->fsource = NULL;
		seq[i]->sorce = seq[i]->esorce = NULL;
		seq[i]->simple_sorce = NULL;
		seq[i]->node_logo = seq[i]->node_sorce = 0;
		seq[i]->in = 1;
		seq[i]->out = 1;
		seq[i]->sub = -1;
		seq[i]->edge_weight = NULL;
		seq[i]->passing_seq = NULL;
		seq[i]->f0 = NULL;
		seq[i]->mismatch_num = 0;
	}
	seq[len_b - 1]->out = 0;
	free(seq[len_b - 1]->next);
	seq[len_b - 1]->next = NULL;

	topo* n1 = n;
	int num1 = n1->len - 1;
	int num2 = len_b - 1;
	int cont = 0;

	int s1 = INT_MIN;
	int s2 = 0;
	int s4 = 0;
	char s5 = 0;
	for (int i = n1->len - 1; i > 0; i--)
	{
		if (n1->sort[i]->out == 0)
		{
			if (s1 <= n1->sort[i]->lastsorce)
			{
				s1 = n1->sort[i]->lastsorce;
				num1 = n1->sort[i]->sub;
			}
			s2++;
		}
		if (s2 >= n->last_node_num)
			break;
	}
	while (num1 != -1 && num2 != -1)
	{
		if (n1->sort[num1]->source[NUM2(num2)] / 42 == 3)
		{
			cont = 0;
			seq[num2]->sorce = (char*)mm_malloc(L * sizeof(char));
			seq[num2]->esorce = (char*)mm_malloc(L * sizeof(char));
			seq[num2]->passing_seq = (char*)malloc(sum * sizeof(char));
			memset(seq[num2]->passing_seq, 0, sum * sizeof(char));
			seq[num2]->passing_seq[num] = 1;
			seq[num2]->edge_weight = (int*)malloc(sizeof(int));
			seq[num2]->edge_weight[0] = 1;
			seq[num2]->f0 = (char*)malloc(sizeof(char));			
			n->unsort[n->len + l] = seq[num2];
			n->unsort[n->len + l]->sub = n->len + l;
			l++;
			if (NUM2(num2 - 1) > 0 && ((n1->sort[num1]->fsource[NUM2(num2)] == 1 || n1->sort[num1]->fsource[NUM2(num2)] == -1) ||((n1->sort[num1]->fsource[NUM2(num2)] == 2 || n1->sort[num1]->fsource[NUM2(num2)] == -2) && n1->sort[num1]->fsource[NUM2(num2 - 1)] < 0)))
			//if (NUM2(num2 - 1) > 0 && n1->sort[num1]->fsource[NUM2(num2)] == 1)
				n1->sort[num1]->source[NUM2(num2 - 1)] = 126;
			num2--;
			continue;
		}
		else if (n1->sort[num1]->source[NUM2(num2)] / 42 == 0)
		{
			cont = 3;
			if (n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub > 0 && ((n1->sort[num1]->esource[NUM2(num2)] <= 42 && n1->sort[num1]->esource[NUM2(num2)] >= -42) || ((n1->sort[num1]->esource[NUM2(num2)] > 42 || n1->sort[num1]->esource[NUM2(num2)] < -42) && n1->sort[n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub]->esource[NUM2(num2)] < 0)))
			//if (n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub > 0 && n1->sort[num1]->esource[NUM2(num2)] <= 42)
			{
				s5 = n1->sort[n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub]->esource[NUM2(num2)] % 42;
				s5 = (s5 >= 0 ? s5 : -s5) - 1;
				n1->sort[n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub]->source[NUM2(num2)] = s5;
			}
			num1 = n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub;
			continue;
		}
		else if (n1->sort[num1]->source[NUM2(num2)] / 42 == 1)
		{
			if (num2 == len_b - 1)
			{
				if (n1->sort[n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub]->source[NUM2(num2 - 1)] / 42 == 1)
				{
					n1->sort[num1]->edge_weight[n1->sort[num1]->source[NUM2(num2)] % 42]++;
					free_node(seq[num2]);
				}
				else
				{
					n1->sort[num1]->in++;
					n1->sort[num1]->pre = (poa**)realloc(n1->sort[num1]->pre, n1->sort[num1]->in * sizeof(poa*));
					n1->sort[num1]->pre[n1->sort[num1]->in - 1] = seq[num2 - 1];

					n1->sort[num1]->edge_weight = (int*)realloc(n1->sort[num1]->edge_weight, n1->sort[num1]->in * sizeof(int));
					n1->sort[num1]->edge_weight[n1->sort[num1]->in - 1] = 1;
					n1->sort[num1]->f0 = (char*)realloc(n1->sort[num1]->f0, n1->sort[num1]->in * sizeof(char));
					seq[num2 - 1]->next[seq[num2 - 1]->out - 1] = n1->sort[num1];
					free_node(seq[num2]);
					seq[num2] = n1->sort[num1];
				}
			}
			else if (num2 == 0)
			{
				free_node(seq[0]->pre[0]);
				if (cont == 1 || cont == 5)
				{
					free_node(seq[num2]);
					seq[num2] = n1->sort[num1];
				}
				else
				{
					n1->sort[num1]->out++;
					n1->sort[num1]->next = (poa**)realloc(n1->sort[num1]->next, n1->sort[num1]->out * sizeof(poa*));
					n1->sort[num1]->next[n1->sort[num1]->out - 1] = seq[num2 + 1];
					seq[num2 + 1]->pre[seq[num2 + 1]->in - 1] = n1->sort[num1];
				}
			}
			else
			{
				if (n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub != -1 && n1->sort[n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub]->source[NUM2(num2 - 1)] / 42 == 1)
				{
					n1->sort[num1]->edge_weight[n1->sort[num1]->source[NUM2(num2)] % 42]++;
				}
				else
				{
					n1->sort[num1]->in++;
					n1->sort[num1]->pre = (poa**)realloc(n1->sort[num1]->pre, n1->sort[num1]->in * sizeof(poa*));
					n1->sort[num1]->pre[n1->sort[num1]->in - 1] = seq[num2 - 1];
					
					n1->sort[num1]->edge_weight = (int*)realloc(n1->sort[num1]->edge_weight, n1->sort[num1]->in * sizeof(int));
					n1->sort[num1]->edge_weight[n1->sort[num1]->in - 1] = 1;
					n1->sort[num1]->f0 = (char*)realloc(n1->sort[num1]->f0, n1->sort[num1]->in * sizeof(char));
					seq[num2 - 1]->next[seq[num2 - 1]->out - 1] = n1->sort[num1];
				}

				if (cont == 1 || cont == 5)
				{
				}
				else
				{
					n1->sort[num1]->out++;
					n1->sort[num1]->next = (poa**)realloc(n1->sort[num1]->next, n1->sort[num1]->out * sizeof(poa*));
					n1->sort[num1]->next[n1->sort[num1]->out - 1] = seq[num2 + 1];
					seq[num2 + 1]->pre[seq[num2 + 1]->in - 1] = n1->sort[num1];////
				}

				free_node(seq[num2]);
				seq[num2] = n1->sort[num1];
			}
			cont = 1;
			n1->sort[num1]->passing_seq[num] = 1;
			num1 = n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub;
			num2--;
		}
		else
		{
			s4 = 0;
			for (int s = 0; s < n1->sort[num1]->mismatch_num; s++)
			{
				if (seq[num2]->base == n1->sort[num1]->mismatch_node[s]->base)
				{
					if (num2 != 0)
					{
						if (n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub != -1)
						{
							if (n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->source[NUM2(num2 - 1)] / 42 == 1)
							{
								for (int ss = 0; ss < n1->sort[num1]->mismatch_node[s]->in; ss++)
									if (n1->sort[num1]->mismatch_node[s]->pre[ss] == n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42])
									{
										n1->sort[num1]->mismatch_node[s]->edge_weight[ss]++;
										s2 = -1;
									}
							}
						}
						if (s2 != -1)
						{
							n1->sort[num1]->mismatch_node[s]->in++;
							n1->sort[num1]->mismatch_node[s]->pre = (poa**)realloc(n1->sort[num1]->mismatch_node[s]->pre, n1->sort[num1]->mismatch_node[s]->in * sizeof(poa*));
							n1->sort[num1]->mismatch_node[s]->pre[n1->sort[num1]->mismatch_node[s]->in - 1] = seq[num2 - 1];

							n1->sort[num1]->mismatch_node[s]->edge_weight = (int*)realloc(n1->sort[num1]->mismatch_node[s]->edge_weight, n1->sort[num1]->mismatch_node[s]->in * sizeof(int));
							n1->sort[num1]->mismatch_node[s]->edge_weight[n1->sort[num1]->mismatch_node[s]->in - 1] = 1;
							seq[num2 - 1]->next[seq[num2 - 1]->out - 1] = n1->sort[num1]->mismatch_node[s];
						}
					}
					s4 = 1;
					if (cont == 1 || cont == 4)
					{
						for (int ss = 0; ss < seq[num2 + 1]->in; ss++)
							if (seq[num2 + 1]->pre[ss] == n1->sort[num1]->mismatch_node[s])
							{
								s4 = 2;
								seq[num2 + 1]->edge_weight[ss]++;
								seq[num2 + 1]->in--;
								seq[num2 + 1]->pre = (poa**)realloc(seq[num2 + 1]->pre, seq[num2 + 1]->in * sizeof(poa*));
								seq[num2 + 1]->edge_weight = (int*)realloc(seq[num2 + 1]->edge_weight, seq[num2 + 1]->in * sizeof(int));
								seq[num2 + 1]->f0 = (char*)realloc(seq[num2 + 1]->f0, seq[num2 + 1]->in * sizeof(char));
							}
					}

					if (s4 == 1 && num2 != len_b - 1)
					{
						seq[num2 + 1]->pre[seq[num2 + 1]->in - 1] = n1->sort[num1]->mismatch_node[s];
						n1->sort[num1]->mismatch_node[s]->out++;
						n1->sort[num1]->mismatch_node[s]->next = (poa**)realloc(n1->sort[num1]->mismatch_node[s]->next, n1->sort[num1]->mismatch_node[s]->out * sizeof(poa*));
						n1->sort[num1]->mismatch_node[s]->next[n1->sort[num1]->mismatch_node[s]->out - 1] = seq[num2 + 1];
					}

					n1->sort[num1]->mismatch_node[s]->passing_seq[num] = 1;
					if (s2 == -1)
						cont = 5;
					else
						cont = 4;
					s2 = 0;
					free_node(seq[num2]);
					seq[num2] = n1->sort[num1]->mismatch_node[s];
				}
			}
			if (s4 == 0)
			{
				cont = 2;
				seq[num2]->sorce = (char*)mm_malloc(L * sizeof(char));
				seq[num2]->esorce = (char*)mm_malloc(L * sizeof(char));
				seq[num2]->passing_seq = (char*)malloc(sum * sizeof(char));//
				memset(seq[num2]->passing_seq, 0, sum * sizeof(char));//
				seq[num2]->passing_seq[num] = 1;//
				seq[num2]->edge_weight = (int*)malloc(sizeof(int));
				seq[num2]->edge_weight[0] = 1;
				seq[num2]->f0 = (char*)malloc(sizeof(char));

				n->unsort[n->len + l] = seq[num2];
				n->unsort[n->len + l]->sub = n->len + l;
				l++;

				n1->sort[num1]->mismatch_num++;
				n1->sort[num1]->mismatch_node[n1->sort[num1]->mismatch_num - 1] = seq[num2];
				seq[num2]->mismatch_num = n1->sort[num1]->mismatch_num;
				seq[num2]->mismatch_node[seq[num2]->mismatch_num - 1] = n1->sort[num1];
				for (int s = 0; s < n1->sort[num1]->mismatch_num - 1; s++)
				{
					n1->sort[num1]->mismatch_node[s]->mismatch_num++;
					n1->sort[num1]->mismatch_node[s]->mismatch_node[n1->sort[num1]->mismatch_num - 1] = seq[num2];
					seq[num2]->mismatch_node[s] = n1->sort[num1]->mismatch_node[s];
				}
			}
			num1 = n1->sort[num1]->pre[n1->sort[num1]->source[NUM2(num2)] % 42]->sub;
			num2--;
		}
	}

	while (num2 > -1)
	{	
		seq[num2]->sorce = (char*)mm_malloc(L * sizeof(char));
		seq[num2]->esorce = (char*)mm_malloc(L * sizeof(char));
		seq[num2]->passing_seq = (char*)malloc(sum * sizeof(char));//
		memset(seq[num2]->passing_seq, 0, sum * sizeof(char));//
		seq[num2]->passing_seq[num] = 1;//		
		seq[num2]->edge_weight = (int*)malloc(sizeof(int));
		seq[num2]->edge_weight[0] = 1;
		seq[num2]->f0 = (char*)malloc(sizeof(char));
		n->unsort[n->len + l] = seq[num2];
		l++;
		num2--;
	}
	n->len = n->len + l;
	n->unsort = (poa**)realloc(n->unsort, n->len * sizeof(poa*));
	n->sort = (poa**)realloc(n->sort, n->len * sizeof(poa*));
	return n1;
}

topo* control(topo* p, char* A, int num, int sum, ThreadPool* pool)
{
	pthread_mutex_init(&mutex, NULL);
	char* seq = readseq(A, p->len);
	unsigned int tsl;
	length1 = strlen(seq);
	length2 = p->len;
	if (p->len % L != 0)
		length2 = p->len + (L - p->len % L);
	tsl = (length1 + length2) / L - 1;

	if (length1 >= length2)
	{
		fmaxtag = length2 / L - 1;
		lmaxtag = length1 / L - 1;
	}
	else
	{
		fmaxtag = length1 / L - 1;
		lmaxtag = length2 / L - 1;
	}

	maxtag = length1 / L - 1;

	for (int i = 0; i < p->len; i++)
	{
		if (p->sort[i]->source)
			mm_free(p->sort[i]->source);
		if (p->sort[i]->esource)
			mm_free(p->sort[i]->esource);
		if (p->sort[i]->fsource)
			mm_free(p->sort[i]->fsource);
		if (p->sort[i]->simple_sorce)
			free(p->sort[i]->simple_sorce);
		p->sort[i]->source = (char*)mm_malloc(length1 * sizeof(char));
		p->sort[i]->esource = (char*)mm_malloc(length1 * sizeof(char));
		p->sort[i]->fsource = (char*)mm_malloc(length1 * sizeof(char));
		p->sort[i]->simple_sorce = (int*)malloc((maxtag + 2) * sizeof(int));
	}	
	s_len = p->len;//temp
	int j = 0;
	for (int i = 0; i < tsl; i++)
	{
		if (i <= fmaxtag)
			j++;
		else if (i <= lmaxtag)
			;
		else
			j--;
		lock = 0;
		for (int l = 0; l < j; l++)
		{
			pack* pa = (pack*)malloc(sizeof(pack));
			pa->i = i;
			pa->j = j;
			pa->l = l;
			pa->p = p;
			pa->num = num;
			pa->seq = seq;
			threadPoolAdd(pool, block_alignment, pa);
			//block_alignment(i, j, l, p);
			//block_alignment(pa);
		}
		while (lock != j) {}
	}
	
	pthread_mutex_destroy(&mutex);
	int last = maxtag + 1;
	node_fuse(p, A, num, sum, last);
	return p;
}
