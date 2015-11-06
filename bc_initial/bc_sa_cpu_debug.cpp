/**
 File name: bc_stl_cpu.cu
 Author: Yuede Ji
 Last update: 11:15 10-24-2015
 Description: STL CPU bc on small graph
    (1) read begin position, csr, weight value from binary file
    (2) betweenness centrality from 0 to others

**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"
#include "graph.h"
#include <queue>

using namespace std;

const char output_file[] = "/home/yuede/small_graph/result_bc/bc_sa_cpu.sa";
const char debug_file[] = "/home/yuede/small_graph/result_bc/bc_sa_cpu.debug";

const int INF = 0x7fffffff;
const int V = 218;

path_t bc[V];
index_t sa[V];
index_t sp_count[V];
path_t dist[V];
void print_result()
{
    FILE * fp = fopen(output_file, "w");
    for(int i=0; i<V; ++i)
    {
        //fprintf(fp, "%d %g\n", i, bc[i]);
        fprintf(fp, "%u\n", sa[i]);
    }
    fclose(fp);
}
void print_debug(graph *g)
{
    FILE * fp = fopen(debug_file, "w");
    for(index_t j=0; j<V; ++j)
    {
        fprintf(fp, "\n%u:\t", j);
        for(index_t i=g->beg_pos[j]; i<g->beg_pos[j+1]; ++i)
        {
            fprintf(fp, "%u %lf\t", g->csr[i], g->weight[i]); 
        }
    }
    fclose(fp);
}
index_t sssp(index_t root, graph *g)
{
    for(index_t i=0; i<g->vert_count; ++i)
        dist[i] = INF;
    dist[root] = 0;
    
    memset(sa, -1, sizeof(index_t)*g->vert_count);
    memset(sp_count, 0, sizeof(index_t)*g->vert_count);
    sp_count[root] = 1;
    sa[root] = 0;
    bool flag = true;
    index_t level = 0;
    while(flag)
    {
        flag = false;
        for(index_t j=0; j<g->vert_count; ++j)
        {
            if(sa[j] != level)
                continue;
            for(index_t i=g->beg_pos[j]; i<g->beg_pos[j+1]; ++i)
            {
                index_t v = g->csr[i];
                //printf("weight[%u] = %g\n", v, g->weight[v]);
                if(dist[v] > dist[j] + g->weight[i])
                {
                    dist[v] = dist[j] + g->weight[i];
                    flag = true;
                    sa[v] = level + 1;
                    sp_count[v] = 0;
                }
                if(dist[v] == dist[j] + g->weight[i])
                    sp_count[v] += sp_count[j];
            }
        }
        ++level;
    }
    return level;
}

void bc_one(index_t root, graph *g, index_t level)
{
    path_t bc_root[g->vert_count];
    memset(bc_root, 0, sizeof(path_t)*V);
    //reverse the order
    for(index_t cur=level-1; cur>=0; --cur)
    {
        //printf("cur = %u\n", cur);
        for(index_t i=0; i<g->vert_count; ++i)
        {
            if(sa[i] == cur)
            {
                index_t w = i;
                //printf("g->beg_pos[%u] = %u\n", w, g->beg_pos[w]);
                //Undirected graph
                for(index_t j=g->beg_pos[w]; j<g->beg_pos[w+1]; ++j)
                {
                    index_t v = g->csr[j];
                    //printf("v = %u\n", v);
                    if(dist[w] == dist[v] + g->weight[j])
                    {
                        if(sp_count[w] != 0)
                            bc_root[v] += sp_count[v]*1.0/sp_count[w]*(1+bc_root[w]); 
                    }
                }
                if(w != root)///cur > 0 can avoid this branch?
                    bc[w] += bc_root[w];
            }
        }
    }
}
void bc_all(graph *g)
{
    
    index_t sigma[g->vert_count];
    
    index_t st[g->vert_count];// stack, record traverse order, from 0 to g->vert_count - 1

    memset(bc, 0, sizeof(path_t)*g->vert_count);

    sssp(0, g);
    /**

    //calculate all the vertex from 0 to g->vert_count, 
    for(index_t i=0; i<g->vert_count; ++i)
    {
        //take each i as source vertex
        index_t depth = sssp(i, g);
        bc_one(i, g, depth);
    }
    **/
    
}
int main(int args, char ** argv)
{
    printf("Input: ./bfs_small_graph /path/to/beg /path/to/csr thread-count\n");

    if(args != 5)
        exit(-1);
    const char *beg_filename = argv[1];
    const char *csr_filename = argv[2];
    const char *weight_filename = argv[3];
    const int thd_count = atoi(argv[4]);
    
    graph *g = new graph(beg_filename, csr_filename, weight_filename);
    
    bc_all(g);// g->vert_count, g->edge_count);
    
    print_result();

    print_debug(g);

    return 0;
}
