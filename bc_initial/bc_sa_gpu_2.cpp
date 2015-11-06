/**
 File name: bc_sa_gpu.cu
 Author: Yuede Ji
 Last update: 17:10 10-26-2015
 Description: GPU bc on small graph
    (1) read begin position, csr, weight value from binary file
    (2) betweenness centrality

**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"
#include "graph.h"

using namespace std;

const char output_file[] = "/home/yuede/small_graph/result_bc/bc_sa_gpu.result";
const char sp_count_file[] = "/home/yuede/small_graph/result_bc/bc_sa_gpu.sp_count";
const char dist_file[] = "/home/yuede/small_graph/result_bc/bc_sa_gpu.dist";
const char sa_file[] = "/home/yuede/small_graph/result_bc/bc_sa_gpu.sa";

const int INF = 0x7fffffff;
const int V = 218;

path_t bc[V];
index_t sa[V];
index_t sp_count[V];
path_t local_dist[V];
void print_result()
{
    FILE * fp = fopen(output_file, "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%d %g\n", i, bc[i]);
    }
    fclose(fp);
}

void print_debug()
{       
     FILE * fp_count = fopen(sp_count_file, "w");
     for(index_t i=0; i<V; ++i)
         fprintf(fp_count, "%u\n", sp_count[i]);
     fclose(fp_count);
 
     FILE * fp_sa = fopen(sa_file, "w");
     for(index_t i=0; i<V; ++i)
         fprintf(fp_sa, "%u\n", sa[i]);
     fclose(fp_sa);
 
     FILE * fp_dist = fopen(dist_file, "w");
     for(index_t i=0; i<V; ++i)
         fprintf(fp_dist, "%g\n", local_dist[i]);
     fclose(fp_dist);
}

/**__global__ void sssp(index_t root, index_t * dev_sssp, index_t *dev_sp_count, index_t * dev_beg_pos, index_t * dev_csr, path_t * dev_weight, path_t * dev_dist, bool *flag)
{
    __shared__ path_t dist[V];
    index_t id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id == 0)
    {
        dist[root] = 0;
        dev_sp_count[root] = 1;
    }
    index_t level = 0;
    __syncthreads();
    
    while(*flag)
    {
        *flag = false;
        for(index_t j=dev_beg_pos[id]; j<dev_beg_pos[id+1]; ++j)
        {
            if(dev_dist[dev_csr[j]] < INF)
            {
                if(dist[id] > dist[dev_csr[j]] + dev_weight[j])
                {
                    dist[id] = dist[dev_csr[j]] + dev_weight[j];
                    if(!(*flag))
                        *flag = true;
                    sa[id] = level + 1;
                    dev_sp_count[id] = 0;
                }
                if(dist[id] == dist[dev_csr[j]] + dev_weight[j])
                    dev_sp_count[id] += dev_sp_count[j];
            }
        }
        ++level;
        __syncthreads();
    }
    __syncthreads();
    dev_dist[id] = dist[id];
    __syncthreads();
}
__global__ void bc_one(index_t root, index_t * dev_beg_pos, index_t * dev_csr, index_t * dev_weight, index_t * dev_sp_count, index_t * dev_bc, index_t * dev_bc_tmp)
{
    index_t id = threadIdx.x + blockIdx.x * blockDim.x;
    
    path_t bc_root[V];
    memset(bc_root, 0, sizeof(path_t)*V);
    //reverse the order
    for(index_t cur=level-1; cur>=0; --cur)
    {
        for(index_t i=0; i<V; ++i)
        {
            if(sa[i] == cur)
            {
                index_t w = i;
                //Undirected graph
                for(index_t j=g->beg_pos[w]; j<g->beg_pos[w+1]; ++j)
                {
                    index_t v = g->csr[j];
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
   
}**/
__global__ void bc_all(index_t * dev_beg_pos, index_t * dev_csr, path_t *dev_weight, path_t * dev_bc, bool *flag)
{
    //index_t id = threadIdx.x + blockIdx.x * blockDim.x;
    
    __shared__ index_t dev_sp_count[V];

    __shared__ path_t dev_bc_tmp[V];

    __shared__ path_t dev_dist[V];

    __shared__ index_t dev_sa[V];

    index_t root = blockIdx.x;
    index_t dest = threadIdx.x;

    dev_dist[dest] = INF;
    dev_sp_count[dest] = 0;
    dev_bc_tmp[dest] = 0;

    if(dest == 0)
    {
        dev_dist[root] = 0;
        dev_sp_count[root] = 1;
    }
    index_t level = 0;
    __syncthreads();
    
    while(*flag)
    {
        *flag = false;
        for(index_t j=dev_beg_pos[dest]; j<dev_beg_pos[dest+1]; ++j)
        {
            if(dev_dist[dev_csr[j]] < INF)
            {
                if(dev_dist[dest] > dev_dist[dev_csr[j]] + dev_weight[j])
                {
                    dev_dist[dest] = dev_dist[dev_csr[j]] + dev_weight[j];
                    if(!(*flag))
                        *flag = true;
                    dev_sa[dest] = level + 1;
                    dev_sp_count[dest] = 0;
                }
                if(dev_dist[dest] == dev_dist[dev_csr[j]] + dev_weight[j])
                    dev_sp_count[dest] += dev_sp_count[j];
            }
        }
        ++level;
        __syncthreads();
    }
    __syncthreads();
    dev_bc[dest] = dev_dist[dest];
    //printf("%g\n", dev_bc[dest]);
    __syncthreads();
}

void bc_gpu_launch(graph *g)
{
    memset(bc, 0, sizeof(path_t)*g->vert_count);
    index_t v = g->vert_count;
    index_t e = g->edge_count;
    cout<<"v = "<<v<<", e = "<<e<<endl;
   
    index_t *dev_beg_pos;
    vertex_t *dev_csr;
    path_t *dev_sa_global;
    path_t *dev_weight;
    path_t *dev_bc;
    cudaMalloc( (void **) &dev_sa_global, v*sizeof(path_t));
    cudaMalloc( (void **) &dev_beg_pos, (v+1)*sizeof(index_t));
    cudaMalloc( (void **) &dev_csr, e*sizeof(vertex_t));
    cudaMalloc( (void **) &dev_weight, e*sizeof(path_t));
    cudaMalloc( (void **) &dev_bc, v*sizeof(path_t));

    cudaMemcpy(dev_beg_pos, g->beg_pos, (v+1)*sizeof(index_t), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_csr, g->csr, e*sizeof(index_t),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_weight, g->weight, e*sizeof(path_t),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_bc, bc, v*sizeof(path_t), cudaMemcpyHostToDevice);
    
    bool * dev_flag_traverse;
    bool flag_traverse = true;
    
    cudaMalloc( (void **) &dev_flag_traverse, 1 * sizeof(bool));
    cudaMemcpy( dev_flag_traverse, &flag_traverse, 1 * sizeof(bool), cudaMemcpyHostToDevice);

    bc_all<<<1, V>>>(dev_beg_pos, dev_csr, dev_weight, dev_bc, dev_flag_traverse);
    
    cudaDeviceSynchronize();
    
    cudaMemcpy(bc, dev_bc, v*sizeof(path_t), cudaMemcpyDeviceToHost);

    for(index_t i=0; i<v; ++i)
        printf("%g\n", bc[i]);


    cudaFree(dev_sa_global);
    cudaFree(dev_beg_pos);
    cudaFree(dev_csr);
    cudaFree(dev_weight);
    cudaFree(dev_bc);
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
    
    bc_gpu_launch(g);// g->vert_count, g->edge_count);
    
    print_result();

    //print_debug();

    return 0;
}
