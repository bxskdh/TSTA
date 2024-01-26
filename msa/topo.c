#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "poa.h"

static inline int tp1(topo* s, poa* p, int subs)
{
	s->sort[subs] = p;
	s->sort[subs]->node_logo = 0;
	s->sort[subs]->sub = subs;

	int c = 0;
	int max = 0;
	int max_i = 0;

	for (int i = 0; i < s->sort[subs]->in; i++)
	{
		if (s->sort[subs]->pre[i]->node_sorce >= 0)
		{
			if (max < s->sort[subs]->edge_weight[i])
			{
				max = s->sort[subs]->edge_weight[i];
				max_i = i;
			}
			else if (max == s->sort[subs]->edge_weight[i] && s->sort[subs]->pre[max_i]->node_sorce <= s->sort[subs]->pre[i]->node_sorce)///
			{
				max = s->sort[subs]->edge_weight[i];
				max_i = i;
			}
		}
	}

	s->sort[subs]->node_sorce = s->sort[subs]->pre[max_i]->node_sorce + max;
	s->sort[subs]->node_base_len = s->sort[subs]->pre[max_i]->node_base_len + 1;
	s->sort[subs]->node_sorce_source = s->sort[subs]->pre[max_i]->sub;

	p->in_temp = -1;

	subs++;

	for (int j = 0; j < p->out; j++)
	{
		p->next[j]->in_temp--;
		if (p->next[j]->in_temp == 0 && p->next[j]->mismatch_num == 0 && p->next[j]->passing != 2)
			subs = tp1(s, p->next[j], subs);
		else if (p->next[j]->in_temp == 0 && p->next[j]->mismatch_num > 0 && p->next[j]->passing != 2)
		{
			for (int s = 0; s < p->next[j]->mismatch_num; s++)
				if (p->next[j]->mismatch_node[s]->in_temp == 0)
					c++;
			if (c == p->next[j]->mismatch_num)
			{
				c = 0;
				subs = tp1(s, p->next[j], subs);
				for (int ss = 0; ss < p->next[j]->mismatch_num; ss++)
					if(p->next[j]->mismatch_node[ss]->in_temp == 0)
						subs = tp1(s, p->next[j]->mismatch_node[ss], subs);
			}
			c = 0;
		}
	}
	return subs;
}

static inline topo* toposort1(topo* s)
{
	int s1 = 0;
	for (int i = 0; i < s->len; i++)
	{
		s->unsort[i]->in_temp = s->unsort[i]->in;
		s->unsort[i]->passing = 0;
		if (s->unsort[i]->out == 0 && s->unsort[i]->mismatch_num > 0)
		{
			for (int j = 0; j < s->unsort[i]->mismatch_num; j++)
				if (s->unsort[i]->mismatch_node[j]->out != 0)
					s1 = 1;
			if (s1 != 1)
				s->unsort[i]->passing = 2;
		}
		s1 = 0;
	}
		
	int c = 0;
	int subs = 0;
	while (subs < s->len)
		for (int i = 0; i < s->len; i++)
		{
			if (s->unsort[i]->in_temp == 0)
			{
				if (s->unsort[i]->mismatch_num == 0)
				{
					subs = tp1(s, s->unsort[i], subs);
					break;
				}
				else if (s->unsort[i]->in_temp == 0 && s->unsort[i]->mismatch_num > 0)
				{
					for (int j = 0; j < s->unsort[i]->mismatch_num; j++)
						if (s->unsort[i]->mismatch_node[j]->in_temp == 0)
							c++;
					if (c == s->unsort[i]->mismatch_num)
					{
						c = 0;
						subs = tp1(s, s->unsort[i], subs);
						for (int ss = 0; ss < s->unsort[i]->mismatch_num; ss++)
							if (s->unsort[i]->mismatch_node[ss]->in_temp == 0)
								subs = tp1(s, s->unsort[i]->mismatch_node[ss], subs);
						break;
					}
					c = 0;
				}
			}

		}

	return s;
}

static inline topo* modify(topo* p)
{
	topo* s;
	int max = INT_MIN;
	int max_i;
	int Max, Max_i;
	Max = Max_i = 0;
	for (int i = 0; i < p->len; i++)
	{
		if (max <= p->sort[i]->node_sorce)
		{
			max = p->sort[i]->node_sorce;
			max_i = p->sort[i]->sub;
		}
	}

	if (p->sort[max_i]->out == 0)
		return p;
	else
	{
		for (int i = 0; i < p->sort[max_i]->out; i++)
		{
			for (int j = 0; j < p->sort[max_i]->next[i]->in; j++)
				if (p->sort[max_i]->next[i]->pre[j]->node_sorce < p->sort[max_i]->node_sorce && p->sort[max_i]->next[i]->pre[j]->node_sorce > 0)
					p->sort[max_i]->next[i]->pre[j]->node_sorce = -p->sort[max_i]->next[i]->pre[j]->node_sorce;
			p->sort[max_i]->next[i]->node_logo = 4;
		}

		for (int i = max_i + 1; i < p->len; i++)
		{
			if (p->sort[i]->node_sorce >= 0 || p->sort[i]->node_logo == 4)
			{
				for (int j = 0; j < p->sort[i]->in; j++)
				{
					if (p->sort[i]->pre[j]->node_sorce >= 0)
					{
						if (Max < p->sort[i]->edge_weight[j])
						{
							Max = p->sort[i]->edge_weight[j];
							Max_i = j;
						}
						else if (Max == p->sort[i]->edge_weight[j] && p->sort[i]->pre[Max_i]->node_sorce <= p->sort[i]->pre[j]->node_sorce)
						{
							Max = p->sort[i]->edge_weight[j];
							Max_i = j;
						}
					}
				}
				p->sort[i]->node_sorce = p->sort[i]->pre[Max_i]->node_sorce + Max;
				p->sort[i]->node_base_len = p->sort[i]->pre[Max_i]->node_base_len + 1;
				p->sort[i]->node_sorce_source = p->sort[i]->pre[Max_i]->sub;
				p->sort[i]->node_logo = 0;
				Max = Max_i = 0;
			}
		}
		s = modify(p);
	}
	return s;
}

static inline int tp(topo* s, poa* p, int subs)
{
	s->sort[subs] = p;
	s->sort[subs]->node_logo = 0;
	s->sort[subs]->sub = subs;
	p->in_temp = -1;

	subs++;
	for (int j = 0; j < p->out; j++)
	{
		if (p->next[j]->out == 0 && p->next[j]->passing == 1 && p->next[j]->in_temp - 1 == 0)
		{
			p->next[j]->in_temp--;
			if (p->next[j]->in_temp == 0)
				subs = tp(s, p->next[j], subs);
		}
	}
	for (int j = 0; j < p->out; j++)
	{
		p->next[j]->in_temp--;
		if (p->next[j]->in_temp == 0 && p->next[j]->passing != 2)
			subs = tp(s, p->next[j], subs);
	}
	return subs;
}

static inline topo* toposort(topo* s)
{
	int s1 = 0;
	for (int i = 0; i < s->len; i++)
	{
		s->unsort[i]->in_temp = s->unsort[i]->in;
		s->unsort[i]->passing = 0;
		if (s->unsort[i]->out == 0 && s->unsort[i]->mismatch_num > 0)
		{
			for (int j = 0; j < s->unsort[i]->mismatch_num; j++)
				if (s->unsort[i]->mismatch_node[j]->out != 0)
				{
					s->unsort[i]->passing = 1;
					s1 = 1;
				}
			if(s1 != 1)
				s->unsort[i]->passing = 2;
		}
		s1 = 0;
	}

	int subs = 0;
	while (subs < s->len)
		for (int i = 0; i < s->len; i++)
			if (s->unsort[i]->in_temp == 0)
			{
				subs = tp(s, s->unsort[i], subs);
				if (subs + s->last_node_num == s->len)
				{
					for (int i = 0; i < s->len; i++)
						if (s->unsort[i]->in_temp == 0)
							subs = tp(s, s->unsort[i], subs);
				}
				break;
			}
	return s;
}

topo* t_sort(topo* g, int num)
{
	topo* s;
	g->last_node_num = 0;
	for (int i = 0; i < g->len; i++)
		if (g->unsort[i]->out == 0)
			g->last_node_num++;
	if (num != 1)
		s = toposort(g);
	else
	{
		s = toposort1(g);
		s = modify(s);
	}
	for (int i = 0; i < s->len; i++)
		s->unsort[i] = s->sort[i];	
	return s;
}
