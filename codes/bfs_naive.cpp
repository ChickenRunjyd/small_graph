/**
 File name: bfs_small_graph.cu
 Author: Yuede Ji
 Last update: 14:10 10-14-2015
 Description: BFS on small graph
    (1) read begin position, csr from binary file
    (2) using shared memory
**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"
#include "graph.h"

using namespace std;

const int INF = 0x7fffffff;

__global__ void traverse_one(int level, int * dev_sa, index_t * dev_beg_pos, index_t * dev_csr)
{
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(dev_sa[id] == level)///node i belongs to current level
    {
    //int j = dev_beg_pos[id];
        for(index_t j=dev_beg_pos[id]; j<dev_beg_pos[id+1]; ++j)
        {
            if(dev_sa[dev_csr[j]] > level + 1)
            {
                dev_sa[dev_csr[j]] = level + 1;
                //printf("%u\n", dev_csr[j]);
                //cout<<dev_csr[j]<<endl;
            }
        }
    }
}

void bfs_sa(int root, graph *g)
{
    index_t v = g->vert_count;
    index_t e = g->edge_count;
    cout<<"v = "<<v<<", e = "<<e<<endl; 
    int sa[v];
    for(int i=0; i<v; ++i)
        sa[i] = INF;
    int level = 0;
    sa[0] = 0;
    bool flag = true; //flag whether current level has nodes
    int *dev_sa;
    index_t *dev_beg_pos;
    vertex_t *dev_csr;
    cudaMalloc( (void **) &dev_sa, v*sizeof(index_t));
    cudaMalloc( (void **) &dev_beg_pos, (v+1)*sizeof(index_t));
    cudaMalloc( (void **) &dev_csr, e*sizeof(vertex_t));
    cudaMemcpy(dev_sa, sa, v*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_beg_pos, g->beg_pos, (v+1)*sizeof(index_t), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_csr, g->csr, e*sizeof(index_t), cudaMemcpyHostToDevice);
    while(flag)
    {
        flag = false;
        traverse_one<<<1, v>>>(level, dev_sa, dev_beg_pos, dev_csr);
        ++level;
        cudaMemcpy(sa, dev_sa, v*sizeof(int), cudaMemcpyDeviceToHost);
        for(int i=0; i<v; ++i)
            if(sa[i] == level)
            {
                flag = true;
                break;
            }
    }
    cudaMemcpy(sa, dev_sa, v*sizeof(int), cudaMemcpyDeviceToHost);

    FILE * fp = fopen("/home/yuede/dataset_small_graph/bfs_sa.result", "w");
    for(int i=0; i<v; ++i)
    {
        fprintf(fp, "%d\n", sa[i]);
    }
    fclose(fp);
    
    cudaFree(dev_sa);
    cudaFree(dev_beg_pos);
    cudaFree(dev_csr);
    
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


    return 0;
}
