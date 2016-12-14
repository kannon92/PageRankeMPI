CXX=mpic++
CXXFLAGS=-O2 -g -std=c++11

OBJS= SparsePower.o
google: page_rank.cc $(OBJS)
	$(CXX) $(CXXFLAGS) -o google page_rank.cc SparsePower.cc
clean:
	rm -f google SparsePower.o page_rank.o
