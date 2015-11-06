/**
 File name: bc_bottom_up_gpu.cu
 Author: Yuede Ji
 Last update: 1:15 10-28-2015
 Description: GPU bc on small graph
    () read begin position, csr, weight value from binary file
    (2) betweenness centrality
    (3) atomic lock

**/

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"
#include "graph.h"

using namespace std;

const char output_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu.bc";
const char sp_count_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu.sp_count";
const char dist_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu.dist";
const char sa_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu.sa";

const int INF = 0x7fffffff;
const int V = 218;
//const int startnum = 2;
//const int blocknum = 3;

path_t bc[V*V];
index_t sa[V];
index_t sp_count[V];
path_t local_dist[V];

index_t sa_global[V*V];
int sp_count_global[V*V];
path_t dist_global[V*V];
path_t bc_global[V*V];

/*void print_result()
{
    FILE * fp = fopen(output_file, "w");
    for(int i=0; i<V*V; ++i)
    {
        fprintf(fp, "%g\n", bc[i]);
        //fprintf(fp, "%d %g\n", i, bc[i]);
    }
    fclose(fp);
}
*/
void print_debug()
{       
    FILE * fp_count = fopen(sp_count_file, "w");
    for(index_t i=0; i<V*V; ++i)
        fprintf(fp_count, "%d\n", sp_count_global[i]);
    fclose(fp_count);
 
    FILE * fp_sa = fopen(sa_file, "w");
    for(index_t i=0; i<V*V; ++i) 
        fprintf(fp_sa, "%u\n", sa_global[i]);
    fclose(fp_sa);
 
    FILE * fp_dist = fopen(dist_file, "w");
    for(index_t i=0; i<V*V; ++i)
        fprintf(fp_dist, "%g\n", dist_global[i]);
    fclose(fp_dist);

    FILE * fp_bc = fopen(output_file, "w");
    for(int i=0; i<V*V; ++i)
    {
        fprintf(fp_bc, "%g\n", bc_global[i]);
        //fprintf(fp, "%d %g\n", i, bc[i]);
    }
    fclose(fp_bc);
}
__global__ void bc_merge(path_t * dev_bc)
{
    int id = threadIdx.x;
    int bc_id = id + blockDim.x;
    while(bc_id < blockDim.x * blockDim.x)
    {
        dev_bc[id] += dev_bc[bc_id];
        bc_id += blockDim.x;
    }
    __syncthreads();
}

__global__ void bc_all(index_t * dev_beg_pos, 
        index_t * dev_csr, 
        path_t *dev_weight, 
        path_t * dev_bc, 
        index_t * dev_sa_global, 
        int * dev_sp_count_global, 
        path_t * dev_dist_global, 
        int start_vert, 
        int end_vert)
{
    //printf("block_dim = %d\n", blockDim.x);
    __shared__ int shared_sp_count[V];
    __shared__ path_t shared_bc_tmp[V];
    __shared__ path_t shared_dist[V];
    __shared__ index_t shared_sa[V];
    __shared__ bool flag;
    __shared__ int level;
    index_t root = blockIdx.x;
    index_t dest = threadIdx.x;
    index_t id = threadIdx.x + blockIdx.x * blockDim.x;
    //initialize
    shared_dist[dest] = INF;
    shared_sp_count[dest] = 0;
    shared_bc_tmp[dest] = 0;
    shared_sa[dest] = INF;
    
    __syncthreads();
    if(root < start_vert || root >= end_vert)
        return;
    
	if(dest == 0)
    {
        shared_dist[root] = 0;
        shared_sp_count[root] = 1;
        level = 0;
        flag = true;
        shared_sa[root] = 0;
        //printf("%u, %u\n", root, id);
    }
    //printf("%u\n", id);
    __syncthreads();
    //Step: sssp 
    //printf("shared_dist[%u] = %g\n", root, shared_dist[root]);
    int iteration=0;
    assert(gridDim.x >= V);
    while(flag)
    {
        __syncthreads();///the problem!
        flag = false;
//        bool flag_one = false;
//        int prev = -1;
//        iteration++;
        __syncthreads();
        for(index_t j=dev_beg_pos[dest]; j<dev_beg_pos[dest+1]; ++j)
        {
            index_t nebr=dev_csr[j];
//            if(shared_dist[dest] > shared_dist[nebr] + dev_weight[j] && shared_sa[nebr] < iteration)
            if(shared_dist[dest] > shared_dist[nebr] + dev_weight[j])
            {
//                shared_sa[dest] = iteration;
                shared_dist[dest] = shared_dist[nebr] + dev_weight[j];
//                shared_sp_count[dest] = 0;
//                shared_sp_count[dest] = shared_sp_count[nebr];
//                prev = nebr;
//                level = iteration;
                flag = true;
//                flag_one = true;
            }
        }
        __syncthreads(); 
    }
    dev_dist_global[id] = shared_dist[dest];
    __syncthreads();
}

//    dev_sa_global[id] = shared_sa[dest];
//    dev_sp_count_global[id] = shared_sp_count[dest];
    /*    if(flag_one)
        {
            for(index_t j=dev_beg_pos[dest]; j<dev_beg_pos[dest+1]; ++j)
            {
                index_t nebr=dev_csr[j];
                if(shared_dist[dest] == shared_dist[dev_csr[j]] + dev_weight[j]) //&& dev_csr[j] != prev)
                {
                    shared_sp_count[dest] += shared_sp_count[dev_csr[j]];
//                     shared_sa[dest] = shared_sa[dev_csr[j]] + ;
                   // if(shared_sa[dest] < shared_sa[dev_csr[j]] + )
                   // {
                   //     shared_sa[dest] = shared_sa[dev_csr[j]] + ;
                   //     if(shared_sa[dest] > level)
                   //         level = shared_sa[dest];
                   // }
                }
            }
        }
        __syncthreads();
    }
//    dev_dist_global[id] = shared_dist[dest];
//    __syncthreads();
    __syncthreads();
    if(shared_sa[dest] < INF -10 && shared_sa[dest] > level)
    {
        //printf("ok\n");
        //level = 0;
        //printf("%u\n", shared_sa[dest]);
        level = shared_sa[dest];
    }
    
    //__syncthreads();
    //if(dest == 0)
    //    printf("level = %d\n", level);
   // dev_bc[id] = shared_sp_count[dest];
    dev_sa_global[id] = shared_sa[dest];
    dev_sp_count_global[id] = shared_sp_count[dest];
    dev_dist_global[id] = shared_dist[dest];

    __syncthreads();
    */
/*
    //printf("\n");
   //Step 2: bc_one
    //printf("%u %d\n", root, level);
    while(level>=0)
    {
        //printf("level = %d\n", level);

        if(shared_sa[dest] == level)
        {
            for(index_t j=dev_beg_pos[dest]; j<dev_beg_pos[dest+1]; ++j)
            {
                if(shared_dist[dev_csr[j]] == shared_dist[dest] + dev_weight[j])
                {
                    ///shared_bc_tmp[dev_csr[j]]  may be changed
                    //printf("level = %d, sa[%u] = %u, sa[%u] = %u\n", level, dest, shared_sa[dest], dev_csr[j], shared_sa[dev_csr[j]]);
                    if(shared_sp_count[dev_csr[j]] != 0)
                        shared_bc_tmp[dest] += shared_sp_count[dest]*(1.0 + shared_bc_tmp[dev_csr[j]])/(shared_sp_count[dev_csr[j]]);
                }
                //__syncthreads();
            }
        }
        __syncthreads();
        if(dest == 0)
        {
            level = level - 1;
            //printf("level = %d\n", level);
        }
        __syncthreads();
    }
    __syncthreads();
    if(dest == root)
        shared_bc_tmp[root] = 0;
    //atomic 4
    dev_bc[id] = shared_bc_tmp[dest];
    
}
*/


void bc_gpu_launch(graph *g, int start_vert, int end_vert)
{
    index_t v = g->vert_count;
    index_t e = g->edge_count;
    //cout<<"v = "<<v<<", e = "<<e<<endl;
    //printf("gpu %d, %d\n", start_vert, end_vert); 
    index_t *dev_beg_pos;
    vertex_t *dev_csr;
    path_t *dev_weight;
    path_t *dev_bc;
    
    index_t *dev_sa_global;
    int *dev_sp_count_global;
    path_t *dev_bc_global;
    path_t *dev_dist_global;


    cudaMalloc( (void **) &dev_beg_pos, (v+1)*sizeof(index_t));
    cudaMalloc( (void **) &dev_csr, e*sizeof(vertex_t));
    cudaMalloc( (void **) &dev_weight, e*sizeof(path_t));
    cudaMalloc( (void **) &dev_bc, V*V*sizeof(path_t));

    cudaMalloc( (void **) &dev_sa_global, V*V*sizeof(index_t));
    cudaMalloc( (void **) &dev_sp_count_global, V*V*sizeof(int));
    cudaMalloc( (void **) &dev_dist_global, V*V*sizeof(path_t));
    cudaMalloc( (void **) &dev_bc_global, V*V*sizeof(path_t));

    
    cudaMemcpy(dev_beg_pos, g->beg_pos, (v+1)*sizeof(index_t), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_csr, g->csr, e*sizeof(index_t),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_weight, g->weight, e*sizeof(path_t),cudaMemcpyHostToDevice);
    cudaMemset(dev_bc, 0, V*V*sizeof(path_t));
    
    cudaMemset(dev_sa_global, 0, V*V*sizeof(index_t));
    cudaMemset(dev_sp_count_global, 0, V*V*sizeof(int));
    cudaMemset(dev_dist_global, 0, V*V*sizeof(path_t));
    cudaMemset(dev_bc_global, 0, V*V*sizeof(index_t));
    
    bc_all<<<V, V>>>(dev_beg_pos, dev_csr, dev_weight, dev_bc_global, dev_sa_global, dev_sp_count_global, dev_dist_global, start_vert, end_vert);//, dev_flag_traverse);
    cudaDeviceSynchronize();
    //bc_merge<<<, V>>>(dev_bc);
    cudaDeviceSynchronize();
    
    cudaMemcpy(sa_global, dev_sa_global, V*V*sizeof(index_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(sp_count_global, dev_sp_count_global, V*V*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(dist_global, dev_dist_global, V*V*sizeof(path_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(bc_global, dev_bc_global, V*V*sizeof(path_t), cudaMemcpyDeviceToHost);

    cudaFree(dev_sa_global);
    cudaFree(dev_beg_pos);
    cudaFree(dev_csr);
    cudaFree(dev_weight);
    cudaFree(dev_bc);
}

int main(int args, char ** argv)
{
    //printf("Input: ./bfs_small_graph /path/to/beg /path/to/csr thread-count\n");

    if(args != 7)
        exit(-1);
    const char *beg_filename = argv[1];
    const char *csr_filename = argv[2];
    const char *weight_filename = argv[3];
    const int thd_count = atoi(argv[4]);
    const int start_vert = atoi(argv[5]);
    const int end_vert = atoi(argv[6]);
    
    graph *g = new graph(beg_filename, csr_filename, weight_filename);
    
    bc_gpu_launch(g, start_vert, end_vert);// g->vert_count, g->edge_count);
    
    //print_result();

    print_debug();

    return 0;
}
