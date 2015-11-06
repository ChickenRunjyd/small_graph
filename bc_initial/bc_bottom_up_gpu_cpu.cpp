/**
 File name: bc_bottom_up_gpu.cu
 Author: Yuede Ji
 Last update: 11:15 10-28-2015
 Description: GPU bc on small graph
    (1) read begin position, csr, weight value from binary file
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

const char output_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu.result";
const char sp_count_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu.sp_count";
const char dist_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu.dist";
const char sa_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu.sa";

const int INF = 0x7fffffff;
const int V = 218;
const int blockdim = 256;
const int blocknum = 1;
path_t bc[blockdim * blockdim];
index_t sa[V];
index_t sp_count[V];
path_t local_dist[V];
void print_result()

const char bc_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu_cpu.bc";
const char sp_count_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu_cpu.sp_count";
const char dist_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu_cpu.dist";
const char sa_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_gpu_cpu.sa";

const int INF = 0x7fffffff;
const int V = 218;

path_t dist[V];//dist[V*V];
index_t sa[V];
index_t sp_count[V];
path_t bc[V*V];//*V];

path_t cpu_dist[V];//dist[V*V];
index_t cpu_sa[V];
index_t cpu_sp_count[V];
path_t cpu_bc[V*V];//*V];

path_t gpu_dist_all[V*V];//dist[V*V];
index_t gpu_sa_all[V*V];
index_t gpu_sp_count_all[V*V];
path_t gpu_bc_all[V*V];

path_t cpu_dist_all[V*V];//dist[V*V];
index_t cpu_sa_all[V*V];
index_t cpu_sp_count_all[V*V];
path_t cpu_bc_all[V*V];
/*
void print_result()
{
    FILE * fp = fopen(output_file, "w");
    for(int i=0; i<V*V; ++i)
    {
        fprintf(fp, "%d %g\n", i, bc[i]);
        //fprintf(fp, "%g\n", dist[i]);
    }
    fclose(fp);
}
*/
void print_debug(graph *g)
{
    FILE * fp_count = fopen(sp_count_file, "w");
    for(index_t i=0; i<V*V; ++i)
    {
        if(cpu_sp_count_all[i] == gpu_sp_count_all[i])
            fprintf(fp_count, "%u %u %u\n", i, cpu_sp_count_all[i],  gpu_sp_count_all[i]);
    }
    fclose(fp_count);

    FILE * fp_sa = fopen(sa_file, "w");
    for(index_t i=0; i<V*V; ++i)
    {
        if(cpu_sa_all[i] == gpu_sa_all[i])
            fprintf(fp_sa, "%u %u %u\n", i, cpu_sa_all[i], gpu_sa_all[i]);
    }
    fclose(fp_sa);

    FILE * fp_dist = fopen(dist_file, "w");
    for(index_t i=0; i<V*V; ++i)
    {
        if(cpu_dist_all[i] == gpu_dist_all[i])
            fprintf(fp_dist, "%u %g %g\n", i, cpu_dist_all[i], gpu_dist_all[i]);
    }
    fclose(fp_dist);

    FILE * fp_bc = fopen(bc_file, "w");
    for(index_t i=0; i<V*V; ++i)
    {
        if(cpu_bc[i] == gpu_bc[i])
            fprintf(fp_bc, "%u %g %g\n", i, cpu_bc[i], gpu_bc[i]);
    }
    fclose(fp_bc);
}

index_t cpu_sssp(index_t root, graph *g)
{
    //cout<<"v = "<<v<<", e = "<<e<<endl;
    memset(cpu_sp_count, 0, sizeof(index_t) * V);
    cpu_sp_count[root] = 1;
    for(index_t i=0; i<V; ++i)
    {
        cpu_dist[i] = INF;
        cpu_sa[i] = INF;
    }
    cpu_dist[root] = 0;
    cpu_sa[root] = 0;

    index_t level = 0;
    bool flag = true;
    while(flag)
    {
        flag = false;    
        for(index_t i=0; i<V; ++i)
        {
            bool flag_one = false;
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(cpu_dist[g->csr[j]] < INF)
                {
                    if(cpu_dist[i] > cpu_dist[g->csr[j]] + g->weight[j])
                    {
                        cpu_dist[i] = cpu_dist[g->csr[j]] + g->weight[j];
                        cpu_sp_count[i] = 0;
                        cpu_sa[i] = cpu_sa[g->csr[j]] + 1;
                        if(sa[i] > level)
                            level = sa[i];
                        if(!flag_one)
                            flag_one = true;
                    }
                }
            }
            if(flag_one)
            {
                flag = true;
                for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
                    if(cpu_dist[i] == cpu_dist[g->csr[j]] + g->weight[j])
                    {
                        cpu_sp_count[i] += cpu_sp_count[g->csr[j]];
                    }
            }

        }
       // printf("level = %u\n", level);
        //printf("sa[%u] = %lf\n", v-1, sa[v-1]);
    }
    for(index_t i=0; i<V; ++i)
    {
        cpu_sa_all[i+root*V] = cpu_sa[i];
        cpu_sp_count_all[i+root*V] = cpu_sp_count[i];
        cpu_dist_all[i+root*V] = cpu_dist[i];
    }
    return level;
}

void bc_cpu_one(index_t root, graph *g, index_t level)
{ 
    path_t bc_tmp[V];
    memset(bc_tmp, 0, sizeof(path_t)*V);
    for(index_t cur=level; cur>=0; --cur)
    {
        for(index_t i=0; i<V; ++i)
        {
            if(cpu_sa[i] == cur)
            {
                for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
                {
                    index_t w = g->csr[j];
                    if(cpu_dist[w] == cpu_dist[i] + g->weight[j])
                    {
                        //if(sp_count[w] != 0)
                        bc_tmp[i] += cpu_sp_count[i]*1.0*(1+bc_tmp[w])/cpu_sp_count[w];
                    }

                }

            }
        }
    }
    bc_tmp[root] = 0;
    //printf("before bc\n");
    for(index_t i=0; i<V; ++i)
    {
        cpu_bc[i+V*root] = bc_tmp[i];
        //dist[i+V*root] = dist
        //bc[i+V*root] = bc_tmp[i];
    }
    
    //printf("%g\n", bc[1]);
}
void bc_cpu_launch(graph *g)
{
    memset(bc, 0, sizeof(path_t)*V);
    //index_t level = sssp(0, g);
    for(index_t i=0; i<blocknum; ++i)//g->vert_count; ++i)
    {
        index_t level = cpu_sssp(i, g);
        //printf("%u, %u\n", i, level);
        bc_cpu_one(i, g, level);
        //printf("%u, %u\n", i, level);
    }
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


__global__ void bc_all(index_t * dev_beg_pos, index_t * dev_csr, path_t *dev_weight, path_t * dev_bc)
{
    //printf("block_dim = %d\n", blockDim.x);
    __shared__ int dev_sp_count[V];//using int for atomic
    __shared__ path_t dev_bc_tmp[V];
    __shared__ path_t dev_dist[V];
    __shared__ index_t dev_sa[V];
    __shared__ bool flag;
    __shared__ int level;
    int cur; 
    index_t root = blockIdx.x;
    index_t dest = threadIdx.x;
    index_t id = threadIdx.x + blockIdx.x * blockDim.x;
    //initialize
    dev_dist[dest] = INF;
    dev_sp_count[dest] = 0;
    dev_bc_tmp[dest] = 0;
    dev_sa[dest] = INF;
    
    __syncthreads();
    
	if(dest == 0)
    {
        dev_dist[root] = 0;
        dev_sp_count[root] = 1;
        level = 0;
        flag = true;
        dev_sa[root] = 0;
        cur = 0;
        //printf("%u, %u\n", root, id);
    }
    //printf("%u\n", id);
    __syncthreads();
    //Step1: sssp 
    //printf("dev_dist[%u] = %g\n", root, dev_dist[root]);
    while(flag)
    {
        //printf("while loop\n");
        flag = false;
        bool flag_one = false;
        //int prev = -1;
        for(index_t j=dev_beg_pos[dest]; j<dev_beg_pos[dest+1]; ++j)
        {
            if(dev_dist[dev_csr[j]] < INF - 10)
            {
                if(dev_dist[dest] > dev_dist[dev_csr[j]] + dev_weight[j])
                {
                    dev_dist[dest] = dev_dist[dev_csr[j]] + dev_weight[j];
                    dev_sp_count[dest] = 0;
                    //dev_sp_count[dest] = dev_sp_count[dev_csr[j]];
                    //dev_sa[dest] = dev_sa[dev_csr[j]] + 1;
                    //prev = dev_csr[j];
                    flag_one = true;
                    //if(dev_sa[dest] > level)
                    //    level = dev_sa[dest];
                }
            }
        }
        __syncthreads();//ensure dev_dist[dev_csr[j]] unchanged
        if(flag_one)
        {
            flag = true;
            for(index_t j=dev_beg_pos[dest]; j<dev_beg_pos[dest+1]; ++j)
            {
                if(dev_dist[dest] == dev_dist[dev_csr[j]] + dev_weight[j])
                {
                    dev_sp_count[dest] += dev_sp_count[dev_csr[j]];
                    if(dev_sa[dest] < dev_sa[dev_csr[j]] + 1)
                    {
                        dev_sa[dest] = dev_sa[dev_csr[j]] + 1;
                        if(dev_sa[dest] > level)
                            level = dev_sa[dest];
                    }
                }
            }
        }
        __syncthreads();
    }
    __syncthreads();

    dev_bc[id] = dev_sp_count[dest];

    __syncthreads();
    //printf("\n");
   //Step 2: bc_one
    //printf("level = %d\n", level);
    /**
    while(level>=0)
    {
        //printf("%d\n", level);
        cur = level;
        if(dev_sa[dest] == cur)
        {
            for(index_t j=dev_beg_pos[dest]; j<dev_beg_pos[dest+1]; ++j)
            {
                if(dev_dist[dev_csr[j]] == dev_dist[dest] + dev_weight[j])
                {
                    ///dev_bc_tmp[dev_csr[j]]  may be changed
                    //printf("level = %d, sa[%u] = %u, sa[%u] = %u\n", level, dest, dev_sa[dest], dev_csr[j], dev_sa[dev_csr[j]]);
                    if(dev_sp_count[dev_csr[j]] != 0)
                        dev_bc_tmp[dest] += dev_sp_count[dest]*(1.0 + dev_bc_tmp[dev_csr[j]])/(dev_sp_count[dev_csr[j]]);
                }
                //__syncthreads();
            }
        }
        level = cur - 1;
        //printf("level = %d\n", cur);
    }
    __syncthreads();
    if(dest == root)
        dev_bc_tmp[root] = 0;
    //atomic 4
    dev_bc[id] = dev_bc_tmp[dest];
    **/
}
void bc_gpu_launch(graph *g)
{
    index_t v = g->vert_count;
    index_t e = g->edge_count;
    cout<<"v = "<<v<<", e = "<<e<<endl;
   
    index_t *dev_beg_pos;
    vertex_t *dev_csr;
    index_t *dev_sa_global;
    path_t *dev_weight;
    path_t *dev_bc;

    index_t *dev_sa_global;
    index_t *dev_sp_count_global;
    index_t *dev_dist_global;
    
    cudaMalloc( (void **) &dev_beg_pos, (v+1)*sizeof(index_t));
    cudaMalloc( (void **) &dev_csr, e*sizeof(vertex_t));
    cudaMalloc( (void **) &dev_weight, e*sizeof(path_t));
    cudaMalloc( (void **) &dev_bc, blockdim*blockdim*sizeof(path_t));
    
    cudaMalloc( (void **) &dev_sa_gloal, blockdim*blockdim*sizeof(index_t));
    cudaMalloc( (void **) &dev_sp_count_global, blockdim*blockdim*sizeof(index_t));
    cudaMalloc( (void **) &dev_dist_global, blockdim*blockdim*sizeof(path_t));

    cudaMemcpy(dev_beg_pos, g->beg_pos, (v+1)*sizeof(index_t), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_csr, g->csr, e*sizeof(index_t),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_weight, g->weight, e*sizeof(path_t),cudaMemcpyHostToDevice);
    cudaMemset(dev_bc, 0, blockdim*blockdim*sizeof(path_t));

    cudaMemset(dev_sa_global, 0, blockdim*blockdim*sizeof(index_t));
    cudaMemset(dev_sp_count_global, 0, blockdim*blockdim*sizeof(index_t));
    cudaMemset(dev_dist_global, 0, blockdim*blockdim*sizeof(index_t));
    //cudaMemcpy(dev_bc, bc, v*sizeof(path_t), cudaMemcpyHostToDevice);
   /** 
    bool * dev_flag_traverse;
    bool flag_traverse = true;
    
    cudaMalloc( (void **) &dev_flag_traverse, 1 * sizeof(bool));
    cudaMemcpy( dev_flag_traverse, &flag_traverse, 1 * sizeof(bool), cudaMemcpyHostToDevice);

    int *global_lock;
    cudaMalloc( (void **) &global_lock, 1 * sizeof(int));
    cudaMemset( global_lock, 0, 1 * sizeof(int));
    ***/
    bc_all<<<blocknum, V>>>(dev_beg_pos, dev_csr, dev_weight, dev_bc, dev_sa_global, dev_sp_count_global, dev_dist_global);//, dev_flag_traverse);
    cudaDeviceSynchronize();
    //bc_merge<<<1, V>>>(dev_bc);
    cudaDeviceSynchronize();
    
    cudaMemcpy(gpu_bc_all, dev_bc, blockdim*blockdim*sizeof(path_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(gpu_sa_all, dev_sa_global, blockdim*blockdim*sizeof(index_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(gpu_sp_count_all, dev_sp_count_global, blockdim*blockdim*sizeof(index_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(gpu_dist_all, dev_dist_global, blockdim*blockdim*sizeof(path_t), cudaMemcpyDeviceToHost);
//    for(index_t i=0; i<v; ++i)
//        printf("%g\n", bc[i]);

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
    const int thread_count = atoi(argv[4]);
    
    graph *g = new graph(beg_filename, csr_filename, weight_filename);
    
    bc_cpu_launch(g);
    
    bc_gpu_launch(g);// g->vert_count, g->edge_count);
    
    //print_result();

    print_debug();

    return 0;
}

