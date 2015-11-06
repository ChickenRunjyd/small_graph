/**
 File name: bc_sa_gpu.cu
 Author: Yuede Ji
 Last update: 17:10 10-26-2015
 Description: GPU bc on small graph
    (1) read begin position, csr, weight value from binary file
    (2) betweenness centrality
    (3) atomic lock

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
/**
struct Lock
{
    int *mutex;
    Lock( void )
    {
        cudaMalloc( (void**)&mutex, sizeof(int) );
        cudaMemset( mutex, 0, sizeof(int) );
    }
 
    ~Lock( void )
    {
        cudaFree( mutex );
    }

    __device__ void lock( void )
    {
        while( atomicCAS( mutex, 0, 1 ) != 0 );
    }
 
    __device__ void unlock( void )
    {
        atomicExch( mutex, 0 );
    }
};
*/
void print_result()
{
    FILE * fp = fopen(output_file, "w");
    for(int i=0; i<V; ++i)
    {
        //fprintf(fp, "%d %g\n", i, bc[i]);
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
/**
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
        __double_as_longlong(val +
        __longlong_as_double(assumed)));

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
        } while (assumed != old);

    return __longlong_as_double(old);
}
**/
__device__ double fatomicMin(double *addr, double value)
{
    double old = *addr, assumed;
    if(old <= value) return old;
    do
    {
        assumed = old;
        old = atomicCAS((unsigned int*)addr, __longlong_as_double(assumed), __longlong_as_double(value));
    }while(old!=assumed);
    return old;
}
__global__ void bc_all(index_t * dev_beg_pos, index_t * dev_csr, path_t *dev_weight, path_t * dev_bc, int * global_lock)
{
    //printf("bc_all\n");
    __shared__ int dev_sp_count[V];//using int for atomic
    __shared__ path_t dev_bc_tmp[V];
    __shared__ path_t dev_dist[V];
    __shared__ index_t dev_sa[V];
    __shared__ bool flag;
    __shared__ index_t level;
    __shared__ index_t cur; 
    __shared__ int block_lock;
    index_t root = blockIdx.x;
    index_t dest = threadIdx.x; 
    //initialize
    dev_dist[dest] = INF;
    dev_sp_count[dest] = 0;
    dev_bc_tmp[dest] = 0;
    dev_sa[dest] = INF;

    if(dest == 0)
    {
        dev_dist[root] = 0;
        dev_sp_count[root] = 1;
        level = 0;
        flag = true;
        dev_sa[root] = 0;
        block_lock = 0;
    }
    __syncthreads();
    //Step1: sssp 
    //printf("dev_dist[%u] = %g\n", root, dev_dist[root]);
    while(flag)
    {
        //printf("while loop\n");
        flag = false;
        cur = level;
        if(dev_sa[dest] == level)
        {
            for(index_t j=dev_beg_pos[dest]; j<dev_beg_pos[dest+1]; ++j)
            {
                if(dev_dist[dev_csr[j]] > dev_dist[dest] + dev_weight[j])
                {
                    dev_sa[dev_csr[j]] = cur + 1;
                    dev_sp_count[dev_csr[j]] = 0;
                    if(!(flag))
                        flag = true;
                    
                    //atomic 1
                    /*while(atomicCAS(&block_lock, 0, 1)!=0)
                    {
                        printf("%d\n", block_lock);
                    }
                    */
                    if(dev_dist[dev_csr[j]] > dev_dist[dest] + dev_weight[j])
                    {
                        dev_dist[dev_csr[j]] = dev_dist[dest] + dev_weight[j];
                    }
                   // printf("%d\n", block_lock);
                   // atomicExch(&block_lock, 0);
                    //printf("%g\n", dev_dist[dev_csr[j]]);
                    /*lock_a.lock();
                    if(dev_dist[dev_csr[j]] > dev_dist[dest] + dev_weight[j])
                    {
                        dev_dist[dev_csr[j]] = dev_dist[dest] + dev_weight[j];
                    }
                    printf("%u\n", dest);
                    */
                    //lock_a.unlock();
                    //dev_dist[dev_csr[j]] = fatomicMin(&dev_dist[dev_csr[j]], dev_dist[dest] + dev_weight[j]);

                }

                if(dev_dist[dev_csr[j]] == dev_dist[dest] + dev_weight[j])
                {
                    //atomic 2
                    dev_sp_count[dev_csr[j]] += dev_sp_count[dest];
                    //atomicAdd(&dev_sp_count[dev_csr[j]], dev_sp_count[dest]);
                }
            }
        }
        level = cur + 1;
        __syncthreads();
    }
    __syncthreads();
    //printf("level = %u\n", level);
    //dev_bc[dest] = dev_dist[dest]; //dev_dist[dest];
    //printf("%g\n", dev_dist[dest]);
    //__syncthreads();
   
    //Step 2: bc_one
    while(level>=0)
    {
        cur = level;
        if(dev_sa[dest] == cur)
        {
            for(index_t j=dev_beg_pos[dest]; j<dev_beg_pos[dest+1]; ++j)
            {
                if(dev_dist[dest] == dev_dist[dev_csr[j]] + dev_weight[j])
                {
                    //atomic 3
                    if(dev_sp_count[dest] != 0)
                        dev_bc_tmp[dev_csr[j]] += dev_sp_count[dest]*(1.0+dev_bc_tmp[dest])/dev_sp_count[dest];
                }
            }
        }
        
        level = cur - 1;
        __syncthreads();
    }

    if(dest == root)
        dev_bc_tmp[root] = 0;
    __syncthreads();
    //atomic 4
    dev_bc[dest] += dev_bc_tmp[dest];
    __syncthreads();
}
/**
__global__ void bc_all(index_t * dev_beg_pos, index_t * dev_csr, path_t *dev_weight, path_t * dev_bc, bool *flag)
{
    index_t root = blockIdx.x;
    index_t dest = threadIdx.x;
    __shared__ path_t s[V];
    s[dest] = dest;
    dev_bc[dest] = s[dest]; 
}
**/
void bc_gpu_launch(graph *g)
{
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
   /** 
    bool * dev_flag_traverse;
    bool flag_traverse = true;
    
    cudaMalloc( (void **) &dev_flag_traverse, 1 * sizeof(bool));
    cudaMemcpy( dev_flag_traverse, &flag_traverse, 1 * sizeof(bool), cudaMemcpyHostToDevice);
**/
    int *global_lock;
    cudaMalloc( (void **) &global_lock, 1 * sizeof(int));
    cudaMemset( global_lock, 0, 1 * sizeof(int));
    bc_all<<<V, V>>>(dev_beg_pos, dev_csr, dev_weight, dev_bc, global_lock);//, dev_flag_traverse);
    
    cudaDeviceSynchronize();
    
    cudaMemcpy(bc, dev_bc, v*sizeof(path_t), cudaMemcpyDeviceToHost);

//    for(index_t i=0; i<v; ++i)
//        printf("%g\n", bc[i]);

    cudaFree(dev_sa_global);
    cudaFree(dev_beg_pos);
    cudaFree(dev_csr);
    cudaFree(dev_weight);
    cudaFree(dev_bc);
    cudaFree(global_lock);
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
