#ifndef PA2_H
#define PA2_H

#include <iostream>
#include <vector>
#include <cmath>
#include <bits/stdc++.h>
#include <stack>
#include <queue>
#include <set>

struct Vertex{


        Vertex(); // Default Constructor
        Vertex(int index, std::vector<double> hobbies); // Constructor
        Vertex & operator=(const Vertex &right); // Overloading assignment operator
        bool operator==(const Vertex &right); // Overloading comparison operator

        int person; // Person in the social graph
        std::vector<double> hobbies; // Vector of hobbies
};

struct Edge{

        Edge(Vertex v1, Vertex v2, double weight); // Constructor

        Vertex v1; // Vertex 1
        Vertex v2; // Vertex 2

        double weight; // Weight of friendship
};

class GraphGenerator{

    public:

        GraphGenerator(); // Default constructor
        GraphGenerator(std::vector<Edge> edgeList); // Constructor
        void addEdge(Edge edge); // Adds an edge to graph
        GraphGenerator & operator=(const GraphGenerator &right); // Overloading Assignment operator
        std::vector<std::pair<Vertex, double>> adjList[101]; // Graph
};

class GraphOperator{

    public:

        GraphOperator(GraphGenerator graph); // Constructor
        double FindAverageDegree(); // Find the average degree
        std::vector<int> FindHighestDegree(); // Find the vertex with the highest degree
        int FindConnectedNumber(std::vector<Edge> edgeList); // Find the number of connected components
        double* DijkstraAlgo(int vertex); // Finds the shortest path
        double eccentricity(int vertex); // Returns the eccentricity of a vertex
        std::vector<std::vector<int>> FindConnectedComponents(std::vector<Edge> edgeList); // Returns a list of the connected components
        std::vector<std::vector<double>> FindConnectedParameters(std::vector<Edge> edgeList); // Find the diameter, radius, and centers of each component.
        std::vector<int> roots(std::vector<Edge> edgeList);
        double FindOpenTriangles(std::vector<Edge> edgeList);
        double FindClosedTriangles(std::vector<Edge> edgeList);
        double FindTrianglesRatio(std::vector<Edge> edgeList); // Find the ratio between the number of open and closed triangles
        int FindClosestNode(); // Find the closest node from x with an interest level of at least t on hobby h.
        int FindHighestInterest(); // Find a person with the highest interest in h.
        std::pair<int, int> FindDistanceRatio(); // Find a pair of nodes x and y whose ratio between hobby distance and graph distance is smallest.

        GraphGenerator graph;

};

#endif