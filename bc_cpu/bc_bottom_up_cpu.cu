/**
 File name: bc_bottom_up_cpu.cu
 Author: Yuede Ji
 Last update: 14:10 10-25-2015
 Description: Bottom up CPU bc on small graph
    (1) read begin position, csr, weight value from binary file
    (2) betweenness centrality in a small graph

**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"
#include "graph.h"

using namespace std;

const char output_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_cpu_debug.bc";
const char sp_count_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_cpu_debug.sp_count";
const char dist_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_cpu_debug.dist";
const char sa_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_cpu_debug.sa";

const int INF = 0x7fffffff;
const int V = 218;

//int startnum = 0;
//int blocknum = 1;

path_t dist[V];//dist[V*V];
index_t sa[V];
index_t sp_count[V];
path_t bc[V*V];//*V];

path_t dist_all[V*V];//dist[V*V];
index_t sa_all[V*V];
index_t sp_count_all[V*V];
/*void print_result()
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
        fprintf(fp_count, "%u\n", sp_count_all[i]);
    fclose(fp_count);

    FILE * fp_sa = fopen(sa_file, "w");
    for(index_t i=0; i<V*V; ++i)
        fprintf(fp_sa, "%u\n", sa_all[i]);
    fclose(fp_sa);

    FILE * fp_dist = fopen(dist_file, "w");
    for(index_t i=0; i<V*V; ++i)
        fprintf(fp_dist, "%g\n", dist_all[i]);
    fclose(fp_dist);
    FILE * fp_bc = fopen(output_file, "w");
    for(int i=0; i<V*V; ++i)
    {
        fprintf(fp_bc, "%g\n", i, bc[i]);
        //fprintf(fp, "%g\n", dist[i]);
    }
    fclose(fp_bc);
}

index_t sssp(index_t root, graph *g)
{
    //cout<<"v = "<<v<<", e = "<<e<<endl;
    memset(sp_count, 0, sizeof(index_t) * V);
    sp_count[root] = 1;
    for(index_t i=0; i<V; ++i)
    {
        dist[i] = INF;
        sa[i] = INF;
    }
    dist[root] = 0;
    sa[root] = 0;

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
                if(dist[g->csr[j]] < INF)
                {
                    if(dist[i] > dist[g->csr[j]] + g->weight[j])
                    {
                        dist[i] = dist[g->csr[j]] + g->weight[j];
                        sp_count[i] = 0;
                        sa[i] = sa[g->csr[j]] + 1;
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
                    if(dist[i] == dist[g->csr[j]] + g->weight[j])
                    {
                        sp_count[i] += sp_count[g->csr[j]];
                    }
            }

        }
        //printf("level = %u\n", level);
        //printf("sa[%u] = %lf\n", v-1, sa[v-1]);
    }
    for(index_t i=0; i<V; ++i)
    {
        sa_all[i+root*V] = sa[i];
        sp_count_all[i+root*V] = sp_count[i];
        dist_all[i+root*V] = dist[i];
    }
    return level;
}

void bc_one(index_t root, graph *g, index_t level)
{ 
    path_t bc_tmp[V];
    memset(bc_tmp, 0, sizeof(path_t)*V);
    for(index_t cur=level; cur>=0; --cur)
    {
        for(index_t i=0; i<V; ++i)
        {
            if(sa[i] == cur)
            {
                for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
                {
                    index_t w = g->csr[j];
                    if(dist[w] == dist[i] + g->weight[j])
                    {
                        //if(sp_count[w] != 0)
                        bc_tmp[i] += sp_count[i]*1.0*(1+bc_tmp[w])/sp_count[w];
                    }

                }

            }
        }
    }
    bc_tmp[root] = 0;
    //printf("before bc\n");
    for(index_t i=0; i<V; ++i)
    {
        bc[i+V*root] = bc_tmp[i];
        //dist[i+V*root] = dist
        //bc[i+V*root] = bc_tmp[i];
    }
    
    //printf("%g\n", bc[1]);
}
void bc_merge()
{
    for(index_t i=0; i<V; ++i)
    {
        index_t id = i + V;
        while(id < V*V)
        {
            bc[i] += bc[id];
            id += V;
        }
    }
}

void bc_all(graph *g, int start_vert, int end_vert)
{
    memset(bc, 0, sizeof(path_t)*V);
    //index_t level = sssp(0, g);
    for(index_t i=start_vert; i<end_vert; ++i)//g->vert_count; ++i)
    {
        if(i>=V)
            break;
        index_t level = sssp(i, g);
        //printf("%u %u\n", i, level);
        bc_one(i, g, level);
        //printf("%u, %u\n", i, level);
    }
    bc_merge();
}
int main(int args, char ** argv)
{
    //printf("Input: ./bfs_small_graph /path/to/beg /path/to/csr thread-count start_vert end_vert\n");

    if(args != 7)
        exit(-1);
    const char *beg_filename = argv[1];
    const char *csr_filename = argv[2];
    const char *weight_filename = argv[3];
    const int thd_count = atoi(argv[4]);
    const int start_vert = atoi(argv[5]);
    const int end_vert = atoi(argv[6]);
    //printf("cpu %d, %d\n", start_vert, end_vert);
    graph *g = new graph(beg_filename, csr_filename, weight_filename);
    
    bc_all(g, start_vert, end_vert);
    //sssp(0, g);// g->vert_count, g->edge_count);
    
    //print_result();

    print_debug(g);

    return 0;
}
