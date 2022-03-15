#include <iostream>
#include <fstream>
#include <string>
#include "PA2.h"

int main(int argc, char **argv){

    std::string line;   // To read each line from code
	int count = 0;    // Variable to keep count of each line
	
	std::ifstream mFile (argv[1]);   
	
	if(mFile.is_open()) 
	{
		while(mFile.peek()!=EOF)
		{
			getline(mFile, line);
			count++;
		}
		mFile.close();
	}

    std::ifstream edges(argv[1]);
    std::ifstream hobbies(argv[2]);

    std::vector<Edge> edgeList;

    for(int i = 0; i < count; i++){
        Vertex v1;
        Vertex v2;
        std::string person1;
        std::string person2;
        std::string weight;

        std::getline(edges, person1, ',');
        std::getline(edges, person2, ',');
        std::getline(edges, weight, '\n');

        v1.person = std::stoi(person1);
        v2.person = std::stoi(person2);

        Edge e(v1, v2, std::stod(weight));

        edgeList.push_back(e);
    }

    GraphGenerator g(edgeList);
    GraphOperator graph(g);

    for(int i = 1; i <= 100; i++){

        if(graph.graph.adjList[i].size() == 0){
            continue;
        }

        for(int j = 0; j < 20; j++){
        std::string hobby;    
            if(j == 19){
                std::getline(hobbies, hobby, '\n');
                graph.graph.adjList[i].at(0).first.hobbies.push_back(std::stod(hobby));
                continue;
            }
            
            std::getline(hobbies, hobby, ',');
            graph.graph.adjList[i].at(0).first.hobbies.push_back(std::stod(hobby));
        }
    }

    std::vector<int> highDeg = graph.FindHighestDegree();
    std::vector<std::vector<double>> parameters = graph.FindConnectedParameters(edgeList);
    std::pair<int, int> distRat = graph.FindDistanceRatio();

    std::cout << "The average degree:" << std::endl;
    std::cout << graph.FindAverageDegree() << std::endl;
    std::cout << "The vertex with the highest degree:" << std::endl;
    for(size_t i = 0; i < highDeg.size(); i++){

        if(i == highDeg.size() - 1){
            std::cout <<highDeg.at(i) << std::endl;
        }else{
            std::cout <<highDeg.at(i) << ", ";
        }
    }
    std::cout << "The number of connected components:" << std::endl;
    std::cout << graph.FindConnectedNumber(edgeList) << std::endl;
    std::cout << "The diameters, radius, and center(s) of each component:" << std::endl;
    for(size_t i = 0; i < parameters.size(); i++){
        for(size_t j = 0; j < parameters.at(i).size(); j++){
            if(j == parameters.at(i).size() - 1){
                std::cout << parameters.at(i).at(j) << std::endl;
            }else{
                std::cout << parameters.at(i).at(j) << ", ";
            }
        }
    }
    std::cout << "The ratio between the number of open and closed triangles:" << std::endl;
    std::cout << graph.FindTrianglesRatio(edgeList) << std::endl;
    std::cout << "The closest node:" << std::endl;
    std::cout << graph.FindClosestNode() << std::endl;
    std::cout << "A closest node with the highest interest:" << std::endl;
    std::cout << graph.FindHighestInterest() << std::endl;
    std::cout << "The pair of nodes x and y:" << std::endl;
    std::cout << distRat.first << ", " << distRat.second << std::endl;

    return 0;
}