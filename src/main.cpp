#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include "inc/graph.h"
#include "inc/common.h"

using namespace std;

int main(){

    Graph graph(0, "datasets/Enron.g");
    std::cout << "The value of V: " << graph.V << std::endl;
    return 0;

}


