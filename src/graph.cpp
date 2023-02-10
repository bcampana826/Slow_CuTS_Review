#include "inc/graph.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include "inc/score.h"

/**
 * This is a constructor method for our Graphs
 * Param mode - whether or not this is the data graph (1) or the query graph (0)
 * Param input_file - the Filename.g graph to read from
 */
Graph::Graph(unsigned int mode, std::string input_file){

    // Initalize the variables
    vector< set<unsigned int> > ns; //neighbor_set
    vector< set<unsigned int> > r_ns; //r_neighbor_set

    // see read_graph_file - V = number of nodes    
    V = read_graph_file(input_file,ns,r_ns);

    // creates the signatures array. Has a number of elements = number of nodes * number of signatures
    // [v1_sig1, v1_sig2, ..., v1_sigN, v2_sig1, ... v2_sigN, ... vM_sigN]
    signatures = new unsigned int[V*Signature_Properties];

    // list of every node in order?
    order_sequence = new unsigned int[V];

    // We have a pragma line here. This is for OpenMP parallization. More Context needed, Commenting out for now
    // #pragma omp parallel for
    
    // for each node, 
    for(int i=0;i<V;++i){
        // save in ( i * 2 + 0 ) and ( i * 2 + 1 ) 
        signatures[i*Signature_Properties+In_degree_offset] = r_ns[i].size(); // iirc these are the sets with the edges - ns[0] = 1 ns[1] = 0 2 3 4 etc etc
        signatures[i*Signature_Properties+Out_degree_offset] = ns[i].size();

    }

    // init two arrays with size V+1
    neighbors_offset = new unsigned int[V+1];
    r_neighbors_offset = new unsigned int[V+1];

    // #pragma omp parallel for

    // init them with 0s first
    for(unsigned int i=0;i<V+1;++i){
        neighbors_offset[i] = 0;
        r_neighbors_offset[i] = 0;
    }

    // oh okay i see how this works now. 
    // start at the beginning, add the size of the of the previous 
    // so [ 0 , ns[0] , ns[0]+ns[1] , ns[0]+ns[1]+ns[2] , etc etc ]
    for(unsigned int i=1;i<V+1;++i){
        neighbors_offset[i] += neighbors_offset[i-1] + ns[i-1].size();
        r_neighbors_offset[i] += r_neighbors_offset[i-1] + r_ns[i-1].size();
    }

    // E = edges - which is equal to the final value of neighbors_offset!!!!!! so cool
    E = neighbors_offset[V];
    neighbors = new unsigned int[neighbors_offset[V]];
    r_neighbors = new unsigned int[r_neighbors_offset[V]];

    // i dont know if this math works out the way theyd like but i hope it does
    AVG_DEGREE = E/V + 2;

    // 
    unsigned int j = 0;
    unsigned int k = 0;

    // so now this is filling the actual neighbor and reverse neighbor arrays
    for(unsigned int i=0;i<V;++i){
        // set of ns[i] which is all the nodes connected to i
        std::set<unsigned int> s = ns[i];

        for(std::set<unsigned int>::iterator p = s.begin();p!=s.end();p++){
            neighbors[j] = *p;
            j++;
        }

        s = r_ns[i];
        for(std::set<unsigned int>::iterator p = s.begin();p!=s.end();p++){
            r_neighbors[k] = *p;
            k++;
        }
    }

    /**
     * for clarity
     * 
     * neighors offset would look something like this = [ 0 , 1 , 5 , 9 ]
     * then neighbors would look like [ 1 , 0 , 2 , 3 , 4 , 1 , 3 , 4 , 5]
     *                                  ^0  ^1              ^2                
     * so crazy but cool
     * 
    */

    // so mode in there case is 1 if data and 0 if query
    // if statement goes off for query graphs only
    if(!mode){
        //sort_search_order(ns,r_ns);
    }

}

// Remember this is only run for query graphs - basically want to get the order of degrees
void Graph::sort_search_order(vector<set<unsigned int>> ns, vector<set<unsigned int>> r_ns)
{

    unsigned int max_out_degree = 0;
    unsigned int idx;

    // These Offsets should work the same way as above in the Graph Constructor
    parents_offset = new unsigned int[V + 1];
    children_offset = new unsigned int[V + 1];

    parents_offset[0] = children_offset[0] = 0;

    // for each vertex, init the offset to 0
    // if the max_out_degree is less than the the max_out_degree of this vertex, save it and set idx to v
    for (unsigned int v = 0; v < V; ++v)
    {
        parents_offset[v + 1] = children_offset[v + 1] = 0;
        if (max_out_degree < signatures[v * Signature_Properties + Out_degree_offset])
        {
            max_out_degree = signatures[v * Signature_Properties + Out_degree_offset];
            idx = v;
        }
    }

    // we want to start with the largest degree vertex
    order_sequence[0] = idx;

    unsigned int inserted_vertexes = 1;

    /***
     * There is an entire .cpp and .h file for this scoring idea,
     * I didn't get to looking through that but we know what this does -
     * gets the order of the query graph
     *
     */
    while (inserted_vertexes < V)
    {
        Score max_score(0, 0, 0);
        for (unsigned int v = 0; v < V; ++v)
        {
            if (binary_search(v, order_sequence, 0, inserted_vertexes))
            {
                continue;
            }
            unsigned int score1 = get_score1(order_sequence, inserted_vertexes, neighbors, neighbors_offset, v);
            unsigned int score2 = get_score2(neighbors, order_sequence, inserted_vertexes, neighbors_offset, v);
            unsigned int score3 = get_score3(order_sequence, inserted_vertexes, neighbors, neighbors_offset, v, V);
            Score temp_score(score1, score2, score3);
            if (compare_score(temp_score, max_score))
            {
                max_score = temp_score;
                idx = v;
            }
        }
        order_sequence[inserted_vertexes++] = idx;
    }

    // Parents and Children sets LIKE IN THE GRAPH CONSTRUCTOR
    vector<set<unsigned int>> P;
    vector<set<unsigned int>> C;
    // init them all with empty sets of ints
    for (unsigned int i = 0; i < V; ++i)
    {
        set<unsigned int> temp1;
        set<unsigned int> temp2;
        P.push_back(temp1);
        C.push_back(temp2);
    }

    // then, for each vertex in order of the order sequence
    for (unsigned int i = 1; i < V; ++i)
    {
        unsigned int v = order_sequence[i];

        // now this loop goes through all the ones we've past already
        // step 1, 0 moves, step 2, 1, step 3, 2 etc etc
        for (unsigned int j = 0; j < i; ++j)
        {

            // get the orevious ones
            unsigned int t_v = order_sequence[j];

            // ns is going out, hence children
            if (ns[v].find(t_v) != ns[v].end())
            {
                C[i].insert(j);
            }
            // reverse is coming in, hence parent
            if (r_ns[v].find(t_v) != r_ns[v].end())
            {
                P[i].insert(j);
            }
        }
    }

    // okay, now that we filled in our array we can make those offsets like before
    // see example in Graph Constructor for how these are structured
    for (unsigned int i = 1; i < V + 1; ++i)
    {
        parents_offset[i] += parents_offset[i - 1] + P[i - 1].size();
        children_offset[i] += children_offset[i - 1] + C[i - 1].size();
    }

    // Now that we know how many parents we need to store, and children, we can make
    // specific arrays for them
    parents = new unsigned int[parents_offset[V]];
    children = new unsigned int[children_offset[V]];
    unsigned int j = 0;
    unsigned int k = 0;

    for (unsigned int i = 0; i < V; ++i)
    {
        std::set<unsigned int> s = P[i];

        // these are the same as above in the Graph constructor
        for (std::set<unsigned int>::iterator p = s.begin(); p != s.end(); p++)
        {
            parents[j] = *p;
            j++;
        }
        s = C[i];
        for (std::set<unsigned int>::iterator p = s.begin(); p != s.end(); p++)
        {
            children[k] = *p;
            k++;
        }
    }
}

/**
 * input file - file to be read
 * neighbors - array of all the left side numbers
 * r_neighbors - array of all the right side numbers??? DIRECTED VS UNDIRECTED GOTTA BE THE ANSWER??
 * 
 * returns number of vertices
*/
unsigned int read_graph_file(string input_file,vector<set<unsigned int> > &neighbors,   vector<set<unsigned int> > &r_neighbors){

    //double load_start = omp_get_wtime();

    ifstream infile;
    infile.open(input_file);
    if(!infile){
        cout<<"load graph file failed "<<endl;
        exit(-1);
    }

    // String to integer - saves in V the first line which is the number of nodes
    unsigned int V  = 0;
    string line;
    getline(infile,line);
    V = stoi(line);

    const std::string delimter = "\t";
    unsigned int line_index = 0;

    // For all verties, make a set for each one on either side and push down
    for(unsigned int i=0;i<V;++i){
        set<unsigned int> temp1;
        set<unsigned int> temp2;
        neighbors.push_back(temp1);
        r_neighbors.push_back(temp2);
    }

    // now were gonna iterate through the file
    while(getline(infile,line)){
        // Find the split point
        auto pos = line.find(delimter);
        if(pos == std::string::npos){
            continue;
        }

        // Split at that point and convert to ints
        int s = stoi(line.substr(0, pos));
        int t = stoi(line.substr(pos + 1, line.size() - pos - 1));

        // then insert these into the contrary set - since they're sets we dont have to worry about over writes
        // for example
        // neighbors[5].insert(6)
        // neighbors[6].insert(5) 
        //
        neighbors[s].insert(t);
        r_neighbors[t].insert(s);
    }

    infile.close();
    //double load_end = omp_get_wtime();
    return V;
}