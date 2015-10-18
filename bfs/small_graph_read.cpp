/**
 File name: small_graph_read.cu
 Author: Yuede Ji
 Last update: 13:30 10-10-2015
 Description: Using status array to implent GPU version of bfs.
    Calculate the shortest distance from 0 to others
**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//Using arrays to implement queue
char fileout[] = "/home/yuede/dataset_small_graph/contact1.result";
char file_v_e[] = "/home/yuede/dataset_small_graph/contact1.v_e";
char file_beg_pos[] = "/home/yuede/dataset_small_graph/contact1_beg.bin";
char file_csr[] = "/home/yuede/dataset_small_graph/contact1_csr.bin";
/**
char filein[] = "/home/yuede/dataset/kron_10_4.dat";// no need
char fileout[] = "/home/yuede/dataset/kron_10_4.gpu.as.result";
char file_v_e[] = "/home/yuede/dataset/kron_10_4.v_e";
char file_beg_pos[] = "/home/yuede/dataset/kron_10_4.beg.pos";
char file_csr[] = "/home/yuede/dataset/kron_10_4.csr";
**/

const int v_num = 256 * 2;
const int e_num = 256*256;
const int INF = 0x7FFFFFFF;
const int threads_num = 256;

int beg_pos[v_num+1];
int csr[e_num];
int sa[v_num];
//load from .dat files, and store in array csr[N*N], beg_pos[N]

int csr_begin();
void bfs_sa(int root, int v, int e);
__global__ void traverse_one(int level, int * dev_sa, int * dev_beg_pos, int * dev_csr)
{
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(dev_sa[id] == level)///node i belongs to current level
    {
        //int j = dev_beg_pos[id];
        for(int j=dev_beg_pos[id]; j<dev_beg_pos[id+1]; ++j)
        {
            if(dev_sa[dev_csr[j]] > level + 1)
            { 
                dev_sa[dev_csr[j]] = level + 1;
                printf("%d\n", dev_csr[j]);
            }
        }                    
    }
}

int main()
{
    csr_begin();

    //bfs_sa(0, v_num, e_num);

    //FILE * fp_out = fopen(fileout, "w");

    //for(int i=0; i<v_num; ++i)
    //     fprintf(fp_out, "%d\n", sa[i]);
    //fclose(fp_out);
    
    return 0;
}
void bfs_sa(int root, int v, int e)
{
    for(int i=0; i<v; ++i)
        sa[i] = INF;
    int level = 0;
    sa[0] = 0;
    bool flag = true; //flag whether current level has nodes
    
    int *dev_sa;
    int *dev_beg_pos;
    int *dev_csr;

    cudaMalloc( (void **) &dev_sa, v*sizeof(int));
    cudaMalloc( (void **) &dev_beg_pos, (v+1)*sizeof(int));
    cudaMalloc( (void **) &dev_csr, e*sizeof(int));

    cudaMemcpy(dev_sa, sa, v*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_beg_pos, beg_pos, (v+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_csr, csr, e*sizeof(int), cudaMemcpyHostToDevice);
    
    while(flag)
    {
        flag = false;
        traverse_one<<<threads_num, threads_num>>>(level, dev_sa, dev_beg_pos, dev_csr);
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

    cudaFree(dev_sa);
    cudaFree(dev_beg_pos);
    cudaFree(dev_csr);

}

int csr_begin()
{
    FILE * fp_beg = fopen(file_beg_pos, "rb");
    fread(&beg_pos, sizeof(beg_pos), 1, fp_beg); 
    fclose(fp_beg);

    int v = v_num;
    for(int j=0; j<v; ++j)
        printf("%d, %d\n", j, beg_pos[j]);
    /*

    i = 0;
    FILE * fp_csr = fopen(file_csr, "rb");
    while(fread(&csr[i], sizeof(int), 1, fp_csr) != EOF)
    {
        ++i;
    }
    fclose(fp_csr);
    */
    //printf("v = %d, e = %d\n", v, e);
}
