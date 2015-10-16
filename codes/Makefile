exe = bfs_small_graph

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
	./bfs_small_graph /home/yuede/dataset_small_graph/contact1_beg.bin /home/yuede/dataset_small_graph/contact1_csr.bin 218

clean:
	rm -rf $(exe) $(objs) 
