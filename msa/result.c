#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "poa.h"

void printf_con(topo* p)
{
	int s = 0;
	FILE* con = fopen("consensus.fa", "w");
	int max = p->sort[p->len - 1]->node_sorce;
	int max_i = p->sort[p->len - 1]->sub;
	for (int i = p->len - 2; i >= 0; i--)
	{
		if (p->sort[i]->out == 0)
		{
			if (max <= p->sort[i]->node_sorce)
			{
				max = p->sort[i]->node_sorce;
				max_i = p->sort[i]->sub;
			}
			s++;
		}
		if (s >= p->last_node_num)
			break;
	}
	int len = p->sort[max_i]->node_base_len;
	char* c = (char*)malloc((len + 1) * sizeof(char));
	c[len] = '\0';
	len--;
	while (len > -1)
	{
		c[len] = p->sort[max_i]->base;
		max_i = p->sort[max_i]->node_sorce_source;
		len--;
	}

	fprintf(con, "consensus:\n");
	fputs(c, con);
	fclose(con);
}

void printf_result(topo* p, int num, FILE* res)
{
	int len = 0;
	char** r = (char**)malloc(num * sizeof(char*));
	for (int i = 0; i < num; i++)
	{
		r[i] = (char*)malloc((p->len + 1) * sizeof(char));
		memset(r[i], '-', p->len);
	}

	for (int i = 0; i < p->len; i++)
		if (p->sort[i]->node_logo == 0)
		{
			for (int j = 0; j < num; j++)
				if (p->sort[i]->passing_seq[j] == 1)
					r[j][len] = p->sort[i]->base;
			for (int l = 0; l < p->sort[i]->mismatch_num; l++)
				for (int j = 0; j < num; j++)
					if (p->sort[i]->mismatch_node[l]->passing_seq[j] == 1)
					{
						r[j][len] = p->sort[i]->mismatch_node[l]->base;
						p->sort[i]->mismatch_node[l]->node_logo = 1;
					}
			len++;
		}
				
	for (int i = 0; i < num; i++)
	{
		r[i][len] = '\0';
		fprintf(res, ">%d\n", i + 1);
		fputs(r[i], res);
		fputs("\n", res);
	}

	printf_con(p);
}
