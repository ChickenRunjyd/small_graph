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

///const int INF = 0x7fffffff;
const int V = 219;

index_t sa[V];

__global__ void traverse_one(int *dev_level, index_t * dev_beg_pos, index_t * dev_csr, index_t * dev_sa_global, bool * dev_flag_traverse)
{
    __shared__ index_t dev_sa[V];

    int id = threadIdx.x + blockIdx.x * blockDim.x;
    /// initialize shared memory to INF, while sa[0] = 0
    if((*dev_level) == 0)
    {
        dev_sa[id] = 0x7fffffff;
        if(id == 0)
            dev_sa[0] = 0;
        //printf("%d\n", 0x7fffffff);
        //printf("dev_sa[%d] = %u\n", id, dev_sa[id]);
        
    }
    else
    {
        dev_sa[id] = dev_sa_global[id];
        //printf("dev_sa[%d] = %u\n", id, dev_sa[id]);
    }
    __syncthreads();
    if(dev_sa[id] == (*dev_level))///node i belongs to current level
    {
        for(index_t j=dev_beg_pos[id]; j<dev_beg_pos[id+1]; ++j)
        {
            /// Do we nedd atomic?
            if(dev_sa[dev_csr[j]] > (* dev_level) + 1)
            {
                dev_sa[dev_csr[j]] = (* dev_level) + 1;
                if(!(*dev_flag_traverse))
                {
                    *dev_flag_traverse = true;
                    //printf("dev_flag = true\n");
                }
            }
        }
    }
    __syncthreads();
    dev_sa_global[id] = dev_sa[id];    
    __syncthreads();
}

void bfs_sa(int root, graph *g)
{
    index_t v = g->vert_count;
    index_t e = g->edge_count;
    cout<<"v = "<<v<<", e = "<<e<<endl;
    
    index_t *dev_beg_pos;
    vertex_t *dev_csr;
    index_t *dev_sa_global;
    cudaMalloc( (void **) &dev_sa_global, v*sizeof(index_t));
    cudaMalloc( (void **) &dev_beg_pos, (v+1)*sizeof(index_t));
    cudaMalloc( (void **) &dev_csr, e*sizeof(vertex_t));
    cudaMemcpy(dev_beg_pos, g->beg_pos, (v+1)*sizeof(index_t), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_csr, g->csr, e*sizeof(index_t), cudaMemcpyHostToDevice);
    
    bool * dev_flag_traverse;
    //bool * dev_flag_end;
    bool * flag_traverse;
    //bool * flag_end;

    int * level;
    int * dev_level;

    cudaHostAlloc((void **) &level, sizeof(int), cudaHostAllocWriteCombined);
    cudaHostGetDevicePointer( &dev_level, level, 0);
    
    cudaHostAlloc((void **) &flag_traverse, sizeof(bool), cudaHostAllocWriteCombined);
    cudaHostGetDevicePointer( &dev_flag_traverse, flag_traverse, 0);
    /**
    cudaHostAlloc((void **) &flag_end, sizeof(bool), cudaHostAllocWriteCombined);
    cudaHostGetDevicePointer( &dev_flag_end, flag_end, 0);
    **/
    * flag_traverse = true;
    * level = 0;
    while(*flag_traverse)
    {
        *flag_traverse = false;
        traverse_one<<<1, V>>>(dev_level, dev_beg_pos, dev_csr, dev_sa_global, dev_flag_traverse);
        cudaDeviceSynchronize();
        ++(*level);
        //if(*flag_traverse)
        //    printf("true\n");
    }
    
    cudaMemcpy(sa, dev_sa_global, v*sizeof(index_t), cudaMemcpyDeviceToHost);

    FILE * fp = fopen("/home/yuede/small_graph/dataset_small_graph/bfs_multi_invoke.result", "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%u\n", sa[i]);
    }
    fclose(fp);
    
    cudaFree(dev_sa_global);
    cudaFree(dev_beg_pos);
    cudaFree(dev_csr);
    cudaFree(dev_flag_traverse);
    cudaFree(dev_level);
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
