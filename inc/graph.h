#include <vector>
#include "common.h"
#include "omp.h"

using namespace std;

class Graph{
public:
    unsigned int V;
    unsigned int E;
    unsigned int AVG_DEGREE = 0;
    unsigned int * neighbors;
    unsigned int * r_neighbors;
    unsigned int * neighbors_offset;
    unsigned int * r_neighbors_offset;
    unsigned int * parents_offset;
    unsigned int * parents;
    unsigned int * children;
    unsigned int * children_offset;
    unsigned int * order_sequence;
    unsigned int * signatures;
    void sort_search_order(vector< set<unsigned int> > ns,vector< set<unsigned int> > r_ns);

    unsigned int read_graph_file(string input_file,vector<set<unsigned int> > &neighbors, vector<set<unsigned int>> &r_neighbors);

    Graph(unsigned int mode,std::string input_file);
    //mode 0 for query graph, mode 1 for data graph
};

