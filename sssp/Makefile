exe = sssp

cc = "$(shell which nvcc)" 
#flags = -I. -fopenmp -march=athlon64 -O3
flags = -I. #-fopenmp -O3
#flags += -std=c++11

ifeq ($(debug), 1)
	flags+= -DDEBUG 
endif

objs = $(patsubst %.cu,%.o,$(wildcard ../../lib/*.cu))\
			$(patsubst %.cu,%.o,$(wildcard *.cu))

deps = $(wildcard ../../lib/*.h) \
				$(wildcard *.h) \
				Makefile

%.o:%.cu $(deps)
	$(cc) -c $< -o $@ $(flags)

$(exe):$(objs)
	$(cc) $(objs) -o $(exe) $(flags)

test:$(exe)
	./sssp /home/yuede/small_graph/dataset_sssp/begin.bin /home/yuede/small_graph/dataset_sssp/adjacent.bin /home/yuede/small_graph/dataset_sssp/weight.bin 218 

test1:$(exe)
	./sssp /home/yuede/small_graph/dataset_sssp/begin.bin /home/yuede/small_graph/dataset_small_graph/adjacent.bin 219

test_back:$(exe)
	./bfs_small_graph /home/yuede/small_graph/dataset_small_graph/cpu_beg_pos_bwd.bin /home/yuede/small_graph/dataset_small_graph/cpu_adj_list.bwd.0.bin 219

test_forward:$(exe)
	./bfs_small_graph /home/yuede/small_graph/dataset_small_graph/cpu_beg_pos_fwd.bin /home/yuede/small_graph/dataset_small_graph/cpu_adj_list.fwd.0.bin 219

clean:
	rm -rf $(exe) $(objs) 
