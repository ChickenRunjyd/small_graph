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
#include <queue>

using namespace std;

const char output_file[] = "/home/yuede/small_graph/dataset_small_graph/bfs_stl_cpu.result";
const char debug_file[] = "/home/yuede/small_graph/dataset_small_graph/bfs_stl_cpu.debug";
const int INF = 0x7fffffff;
const int V = 219;

queue<index_t> q;

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
void bfs(index_t root, graph *g)
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
            if(sa[v]!=INF)
                continue;
            sa[v] = sa[bottom] + 1;
            q.push(v);
        }
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
    
    bfs(0, g);

    print_result();

    print_debug(g);

    return 0;
}
