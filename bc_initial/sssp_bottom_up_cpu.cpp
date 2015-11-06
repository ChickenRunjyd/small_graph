/**
 File name: sssp_bottom_up_cpu.cu
 Author: Yuede Ji
 Last update: 20:45 10-17-2015
 Description: Bottom up CPU sssp on small graph
    (1) read begin position, csr, weight value from binary file
    (2) single source shortest path from 0 to others

**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"
#include "graph.h"

using namespace std;

const char output_file[] = "/home/yuede/small_graph/dataset_sssp/sssp_bottom_up_cpu.result";
const char debug_file[] = "/home/yuede/small_graph/dataset_sssp/sssp_bottom_up_cpu.debug";

const int INF = 0x7fffffff;
const int V = 218;

path_t sa[V];

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
        fprintf(fp, "\n%g:\t", j);
        for(index_t i=g->beg_pos[j]; i<g->beg_pos[j+1]; ++i)
        {
            fprintf(fp, "%u %g\t", g->csr[i], g->weight[i]); 
        }
    }
    fclose(fp);
}
void sssp(int root, graph *g)
{
    index_t v = g->vert_count;
    index_t e = g->edge_count;

    cout<<"v = "<<v<<", e = "<<e<<endl;

    sa[root] = 0;
    for(index_t i=1; i<v; ++i)
    {
        sa[i] = INF;
    }
    
    int level = 0;
    bool flag = true;
    while(flag)
    {
        flag = false;    
        for(index_t i=1; i<v; ++i)
        {
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(sa[g->csr[j]] < INF)
                {
                    if(sa[i] > sa[g->csr[j]] + g->weight[j])
                    {
                        sa[i] = sa[g->csr[j]] + g->weight[j];
                        if(!flag)
                            flag = true;
                    }
                }
            }
        }
        //printf("level = %d\n", level++);
        //printf("sa[%u] = %lf\n", v-1, sa[v-1]);
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
