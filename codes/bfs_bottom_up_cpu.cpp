/**
 File name: bfs_bottom_up_cpu.cu
 Author: Yuede Ji
 Last update: 14:46 10-16-2015
 Description: Bottom up CPU BFS on small graph
    (1) read begin position, csr from binary file
**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"
#include "graph.h"

using namespace std;

const char output_file[] = "/home/yuede/small_graph/dataset_small_graph/bfs_bottom_up_cpu.result";
const char debug_file[] = "/home/yuede/small_graph/dataset_small_graph/bfs_bottom_up_cpu.debug";
const int INF = 0x7fffffff;
const int V = 219;


index_t sa[V];

index_t imin(index_t a, index_t b)
{
    return a<b?a:b;
}

void print_result()
{
    FILE * fp = fopen(output_file, "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%u\n", sa[i]);
    }
    fclose(fp);
}
void print_debug(graph *g)
{
    FILE * fp = fopen(debug_file, "w");
    for(index_t j=0; j<=16; ++j)
    {
        fprintf(fp, "\n%u:\t", j);
        for(index_t i=g->beg_pos[j]; i<g->beg_pos[j+1]; ++i)
        {
            fprintf(fp, "%u\t", g->csr[i]); 
        }
    }
    fclose(fp);
}
void bfs_sa(int root, graph *g)
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
                //denotes g->csr[j] is visited
                if(sa[g->csr[j]] != INF)
                {
                    if(sa[i] > sa[g->csr[j]] + 1)
                    {
                        sa[i] = sa[g->csr[j]] + 1;
                        if(!flag)
                            flag = true;
                    }
                }
            }

        }
        printf("level = %d\n", level++);
        printf("sa[%u] = %u\n", v-1, sa[v-1]);
    }
}

int main(int args, char ** argv)
{
    printf("Input: ./bfs_small_graph /path/to/beg /path/to/csr thread-count\n");

    if(args != 4)
        exit(-1);
    const char *beg_filename = argv[1];
    const char *csr_filename = argv[2];
    const int thd_cound = atoi(argv[3]);
    
    graph *g = new graph(beg_filename, csr_filename);
    
    bfs_sa(0, g);// g->vert_count, g->edge_count);
    
    print_result();

    print_debug(g);

    return 0;
}
