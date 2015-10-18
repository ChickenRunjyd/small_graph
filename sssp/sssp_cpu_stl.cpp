/**
 File name: sssp_cpu_stl.cu
 Author: Yuede Ji
 Last update: 21:37 10-17-2015
 Description: STL CPU sssp on small graph
    (1) read begin position, csr, weight value from binary file
    (2) single source shortest path from 0 to others

**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"
#include "graph.h"
#include <queue>

using namespace std;

const char output_file[] = "/home/yuede/small_graph/dataset_sssp/sssp_cpu_stl.result";
const char debug_file[] = "/home/yuede/small_graph/dataset_sssp/sssp_cpu_stl.debug";

const int INF = 0x7fffffff;
const int V = 218;

path_t sa[V];
queue<index_t> q;

void print_result()
{
    FILE * fp = fopen(output_file, "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%g\n", sa[i]);
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
void sssp(index_t root, graph *g)
{
    for(index_t i=0; i<g->vert_count; ++i)
        sa[i] = INF;
    q.push(root);
    sa[root] = 0;
    
    while(!q.empty())
    {
        index_t bottom = q.front();
        q.pop();
        for(index_t i=g->beg_pos[bottom]; i<g->beg_pos[bottom+1]; ++i)
        {
            index_t v = g->csr[i];
            //printf("weight[%u] = %g\n", v, g->weight[v]);
            if(sa[v] > sa[bottom] + g->weight[i])
            {
                sa[v] = sa[bottom] + g->weight[i];
                if(v == 1)
                    printf("sa[%u] + %g = %g\n", bottom, g->weight[i], sa[v]);
                if(v == 16)
                    printf("sa[%u] + %g = %g\n", bottom, g->weight[i], sa[v]); 
                q.push(v);
            }
        }
        //printf("level\n");
     }
}
int main(int args, char ** argv)
{
    printf("Input: ./bfs_small_graph /path/to/beg /path/to/csr thread-count\n");

    if(args != 5)
        exit(-1);
    const char *beg_filename = argv[1];
    const char *csr_filename = argv[2];
    const char *weight_filename = argv[3];
    const int thd_cound = atoi(argv[4]);
    
    graph *g = new graph(beg_filename, csr_filename, weight_filename);
    
    sssp(0, g);// g->vert_count, g->edge_count);
    
    print_result();

    print_debug(g);

    return 0;
}
