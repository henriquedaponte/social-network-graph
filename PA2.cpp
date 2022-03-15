#include "PA2.h"


Vertex::Vertex(){
    this->person = 0;
}

Vertex::Vertex(int index, std::vector<double> hobbies){

    this->person = index;

    for(size_t i = 0; i < hobbies.size(); i++){
        
        (this->hobbies).push_back(hobbies.at(i));
    }
}

Vertex & Vertex::operator=(const Vertex &right){

    if (&right == this) 
        return (*this);

    this->person = (right.person);
    
    for(size_t i = 0; i < (right.hobbies).size(); i++){

        if(!(this->hobbies).at(i)){
            (this->hobbies).push_back((right.hobbies).at(i));
        }else{
            (this->hobbies).at(i) = (right.hobbies).at(i);
        }
    }

    return (*this); 
}

bool Vertex::operator==(const Vertex &right){

    bool isEqual = true;

    if(this->person != right.person){
        isEqual = false;
    }

    if((this->hobbies).size() != (right.hobbies).size()){
        isEqual = false;
    }else{

        for(size_t i = 0; i < (this->hobbies).size(); i++){

            if((this->hobbies).at(i) != (right.hobbies).at(i)){
                isEqual = false;
            }
        }
    }

    return isEqual;
}

Edge::Edge(Vertex v1, Vertex v2, double weight){

    this->v1 = v1;
    this->v2 = v2;
    this->weight = weight;
}

GraphGenerator::GraphGenerator(){

    Vertex null;
    for(int i = 0; i <= 100; i++){
        (((this->adjList)[i])).push_back({null, 0});
    }
}

void GraphGenerator::addEdge(Edge edge){

    // In case node is an island
    if(edge.v1 == edge.v2){

        if((adjList[(edge.v1).person]).size() == 0){
            (adjList[(edge.v1).person]).push_back({edge.v1, 0});
        }
    }
    else{

        // Creating head node for v1
        if((adjList[(edge.v1).person]).size() == 0){
            (adjList[(edge.v1).person]).push_back({edge.v1, 0});
        }

        // Creating head node for v2
        if((adjList[(edge.v2).person]).size() == 0){
            (adjList[(edge.v2).person]).push_back({edge.v2, 0});
        }

        (adjList[(edge.v1).person]).push_back({edge.v2, edge.weight}); // Adding v2 to v1's list
        (adjList[(edge.v2).person]).push_back({edge.v1, edge.weight}); // Adding v1 to v2's list
    }
}

GraphGenerator::GraphGenerator(std::vector<Edge> edgeList){

    // Adding all edges into adjency list
    for(size_t i = 0; i < edgeList.size(); i++){
        this->addEdge(edgeList.at(i));
    }
}

GraphGenerator & GraphGenerator::operator=(const GraphGenerator &right){

    if (&right == this) 
        return (*this);

    for(int i = 0; i <= 100; i++){
        this->adjList[i] = right.adjList[i];
    }

    return (*this);
}

GraphOperator::GraphOperator(GraphGenerator graph){

    this->graph = graph;
}

double GraphOperator::FindAverageDegree(){  

    Vertex null;
    double degreeSum = 0.0;
    double nodes = 0.0;

    for(int i = 1; i <= 100; i++){

        degreeSum += (((this->graph).adjList[i]).size() - 1); // Edge with itself does not count
        nodes ++;
    }

    // Rounding up to two decimal digits
    double avg = degreeSum / nodes;
    avg = std::ceil(avg * 100.0) / 100.0;

    return avg;
}

std::vector<int> GraphOperator::FindHighestDegree(){

    int highDeg = 0;
    std::vector<int> highPer;
    highPer.push_back(0);

    for(int i = 1; i <= 100; i++){

        int currDeg = ((this->graph).adjList[i]).size();

        if(currDeg == highDeg){
            highPer.push_back(i);
        }

        if(currDeg > highDeg){
            highDeg = currDeg;
            highPer.at(0) = i;

            if(highPer.size() > 1){
                int j = highPer.size();
                while(j > 1){
                    highPer.pop_back();
                    j--;
                }
            }
        }
    }

    return highPer;
}

// FindConnectedNumber helper
int Find(int n1, int par[], int rank[]){

    // parent 1 1 3 4 5
    //n1      1 2 3 4 5

    int res = n1;
     while(res != par[res]){
        par[res] = par[par[res]];
        res = par[res];
     }

    return res;
}

// FindConnectedNumber helper
int Union(int n1, int n2, int par[], int rank[]){

    int p1 = Find(n1, par, rank);
    int p2 = Find(n2, par, rank);

    if(p1 == p2){
        return 0;
    }

    if(rank[p2] > rank[p1]){
        par[p1] = p2;
        rank[p2] += rank[p1];
    }else{
        par[p2] = p1;
        rank[p1] += rank[p2];
    }

    return 1;
}

int GraphOperator::FindConnectedNumber(std::vector<Edge> edgeList){


    int connectComp = 100;

    // Initializing parent array
    int parent[101];
    for(int i = 1; i <= 100; i++){
        parent[i] = i;
    }

    // Initializing rank array
    int rank[101];
    for(int i = 1; i <= 100; i++){
        rank[i] = 1;
    }
    
    
    for(size_t i = 0; i < edgeList.size(); i++){

        int n1 = ((edgeList.at(i)).v1).person;
        int n2 = ((edgeList.at(i)).v2).person;

        connectComp -= Union(n1, n2, parent, rank);
    }

    return connectComp;
}

double* GraphOperator::DijkstraAlgo(int vertex){


    static double dist[101]; // array to calculate the minimum distance for each node                             
    
    // Initializing arays
    for(int i = 1; i <= 100; i++){

        dist[i] = INT_MAX;   
    }

    std::priority_queue<int> pq; // BFS

    pq.push(vertex);
    dist[vertex] = 0;

    while(!pq.empty()){
        
        int u = pq.top();
        pq.pop();

        for(size_t i = 1; i < (this->graph).adjList[u].size(); i++){

            int v = (this->graph).adjList[u].at(i).first.person;
            double weight = (this->graph).adjList[u].at(i).second;

            if(dist[v] > dist[u] + weight){
                dist[v] = dist[u] + weight;
                pq.push(v);
            }
            
        }
    }
    return dist;
}

double GraphOperator::eccentricity(int vertex){

    double maxDist = 0;
    double* Dij = DijkstraAlgo(vertex);

    for(int i = 1; i <= 100; i++){

        if(*(Dij + i) > maxDist && *(Dij + i) != INT_MAX){

            maxDist = *(Dij + i);
        }
    }

    return maxDist;
}

std::vector<std::vector<int>> GraphOperator::FindConnectedComponents(std::vector<Edge> edgeList){

    std::vector<std::vector<int>> connectedComps;
    std::vector<int> notVisited;

    for(int i = 1; i <= 100; i++){
        notVisited.push_back(i);
    }

    int it = 1;

    while(!notVisited.empty()){

        int vertex = notVisited.at(0);
        std::vector<int> component;
        component.push_back(vertex);
        
        notVisited.erase(std::find(notVisited.begin(), notVisited.end(), vertex));

        std::priority_queue<int> pq; // BFS

        pq.push(vertex);

        while(!pq.empty()){
                
            int u = pq.top();
            pq.pop();

            for(size_t j = 1; j < (this->graph).adjList[u].size(); j++){

                int v = (this->graph).adjList[u].at(j).first.person;

                if(std::find(component.begin(), component.end(), v) == component.end()){
                    
                    component.push_back(v);
                    pq.push(v);
                    notVisited.erase(std::find(notVisited.begin(), notVisited.end(), v));
                }
            }
        }

        connectedComps.push_back(component);
        it++;
    }

    return connectedComps;
}

bool sortVec(const std::vector<double> &a, const std::vector<double> &b){

    if(a.at(0) == b.at(0)){
        return (a.at(1) < b.at(1));
    }

    return (a.at(0) < b.at(0));
    
}

std::vector<std::vector<double>> GraphOperator::FindConnectedParameters(std::vector<Edge> edgeList){

    std::vector<std::vector<double>> parameters;
    std::vector<int> centers;
    double radius;
    double diameter;
    std::vector<double> null;
    std::vector<std::vector<int>> connectedComps = FindConnectedComponents(edgeList);

    for(size_t i = 0; i < connectedComps.size(); i++){

        double maxEcc = -1;
        double minEcc = INT_MAX;


        // Initializing parameters list for each component
        parameters.push_back(null);

        // Computing Diameter
        for(size_t j = 0; j < connectedComps.at(i).size(); j++){

            if(eccentricity(connectedComps.at(i).at(j)) > maxEcc){
                maxEcc = eccentricity(connectedComps.at(i).at(j));
            }
        }

        diameter = maxEcc;
        parameters.at(i).push_back(diameter);

        // Computing radius
        for(size_t j = 0; j < connectedComps.at(i).size(); j++){

            if(eccentricity(connectedComps.at(i).at(j)) == minEcc){
                centers.push_back(connectedComps.at(i).at(j));
            }

            if(eccentricity(connectedComps.at(i).at(j)) < minEcc){
                minEcc = eccentricity(connectedComps.at(i).at(j));
                
                size_t k = centers.size();

                while(k > 0){
                    centers.pop_back();
                    k--;
                }

                centers.push_back(connectedComps.at(i).at(j));
            }
        }

        radius = minEcc;
        parameters.at(i).push_back(radius);

        for(size_t j = 0; j < centers.size(); j++){
            parameters.at(i).push_back(centers.at(j));
        }
    }

    // Sorting
    std::sort(parameters.begin(), parameters.end(), sortVec);

    return parameters;
}

std::vector<int> GraphOperator::roots(std::vector<Edge> edgeList){

    int connectComp = 100;
    std::vector<std::vector<int>> connectedComps;

    // Initializing parent array
    int parent[101];
    for(int i = 1; i <= 100; i++){
        parent[i] = i;
    }

    // Initializing rank array
    int rank[101];
    for(int i = 1; i <= 100; i++){
        rank[i] = 1;
    }
    
    for(size_t i = 0; i < edgeList.size(); i++){

        int n1 = ((edgeList.at(i)).v1).person;
        int n2 = ((edgeList.at(i)).v2).person;

        connectComp -= Union(n1, n2, parent, rank);
    }

    std::vector<int> roots;
    
    bool added = false;
    
    for(int i = 1; i <= 100; i++){

        if(added == false){
            roots.push_back(parent[i]);
        }

        added = false;

        for(size_t j = 0; j < roots.size(); j++){
            if(parent[i] == roots.at(j)){
                added = true;
                break;
            }
        }
    }

    return roots;
}

bool sorter(const Edge &a, const Edge &b){

    if(a.v1.person == b.v1.person){
        return (a.v2.person < b.v2.person);
    }

    return (a.v1.person < b.v1.person);
}

double GraphOperator::FindOpenTriangles(std::vector<Edge> edgeList){

    std::vector<std::unordered_set<int>> openTri;

    for(int i = 1; i <= 100; i++){
        for(size_t j = 1; j < (this->graph).adjList[i].size(); j++){

            int n =  (this->graph).adjList[i].at(j).first.person;

            for(size_t k = 1; k < (this->graph).adjList[n].size(); k++){

                int v = (this->graph).adjList[n].at(k).first.person;

                if(v == i){

                    break;
                }else{
                    
                    std::unordered_set<int> temp;
                    temp.insert(i);
                    temp.insert(n);
                    temp.insert(v);

                    bool added = false;

                    for(size_t h = 0; h < openTri.size(); h++){
                        if(temp == openTri.at(h)){
                            added = true;
                        }
                    }

                    if(!added){
                        openTri.push_back(temp);
                    }
                }
            }
        }
    }

    return openTri.size();
}

double GraphOperator::FindClosedTriangles(std::vector<Edge> edgeList){

    double closedTri = 0.0;

    // Swap
    for(size_t i = 0; i < edgeList.size(); i++){

        if(edgeList.at(i).v1.person > edgeList.at(i).v2.person){

            std::swap(edgeList.at(i).v1.person, edgeList.at(i).v2.person);
        }
    }

    // Sort
    std::sort(edgeList.begin(), edgeList.end(), sorter);

    // Find closed triangles
    for(size_t i = 0; i < edgeList.size() - 1; i++){
        for(size_t j = i + 1; j < edgeList.size(); j++){
            Edge e1 = edgeList.at(i);
            Edge e2 = edgeList.at(j);

            if(e1.v1.person == e2.v1.person){
                for(size_t k = 1; k < (this->graph).adjList[e1.v2.person].size(); k++){
                    if((this->graph).adjList[e1.v2.person].at(k).first.person == e2.v2.person){
                        closedTri += 1.0;
                    }
                }
            }
        }
    }

    return closedTri;
}

double GraphOperator::FindTrianglesRatio(std::vector<Edge> edgeList){

    double openTri = FindOpenTriangles(edgeList);
    double closedTri = FindClosedTriangles(edgeList);
    openTri -= closedTri;

    // Rounding up to four decimal digits
    double ratio = openTri / closedTri;
    ratio = std::ceil(ratio * 10000.0) / 10000.0;

    return ratio;
}

int GraphOperator::FindClosestNode(){

    int x = 39; // Designated node
    int h = 6; // Designated hobby(7) index
    double t = 0.5; // Designated Interest level
    double* Dij = DijkstraAlgo(x); // Minimum distances from designated node
    std::vector<std::pair<int, double>> aboveT; // List of nodes with interest above t on hobby at h and their minimum distance from x
    std::vector<int> minAboveT;

    // Making list aboveT
    for(int i = 1; i <= 100; i++){

        if((this->graph).adjList[i].at(0).first.hobbies.at(h) >= t){
            aboveT.push_back({i,*(Dij + i)});
        }
    }

    double minDist = INT_MAX; // Shortext distance from x
    int idx; // Index of node with minDist and interest level on hobby at h above t

    // Find the closest node to x with interest above t on hobby at h
    for(size_t i = 0; i < aboveT.size(); i++){
        if(aboveT.at(i).second < minDist){
            minDist = aboveT.at(i).second;
            idx = aboveT.at(i).first;
        }
    }

    return idx; // Return idx of closest node
}

int GraphOperator::FindHighestInterest(){

    int h = 6; // Index of hobby 7
    double highInt = 0; // highest interest on h quantified
    int idx; // Index of person with highest interest on h

    // Find index of person with highest interest on h
    for(int i = 1; i <= 100; i++){
        if((this->graph).adjList[i].at(0).first.hobbies.at(h) > highInt){

            highInt = (this->graph).adjList[i].at(0).first.hobbies.at(h);
            idx = i;
        }
    }

    // Returns index of person with highest interest on h
    return idx;
}

double hobbyDistance(Vertex v1, Vertex v2){

    double hobbyDist = 0;

    for(size_t i = 0; i < v1.hobbies.size(); i++){

        double dif = v1.hobbies.at(i) - v2.hobbies.at(i);
        hobbyDist += pow(dif, 2);
    }

    hobbyDist = sqrt(hobbyDist);

    return hobbyDist;
}

std::pair<int, int> GraphOperator::FindDistanceRatio(){

    int idx1;
    int idx2;

    std::pair<int, int> idxs;

    double distHobbyRatio = INT_MAX;


    for(int i = 1; i <= 100; i++){

        double* Dij = DijkstraAlgo(i);

        for(int j = 1; j <= 100; j++){

            double hobbyDist = hobbyDistance((this->graph).adjList[i].at(0).first, (this->graph).adjList[j].at(0).first);
            double dist = *(Dij + j);

            // Adressing edge case
            if(dist == INT_MAX || dist == 0){
                continue;
            }

            double ratio = hobbyDist / dist;

            if(ratio < distHobbyRatio){
                distHobbyRatio = ratio;
                idx1 = i;
                idx2 = j;
            }

        }
    }

    idxs = {idx1, idx2};

    if(idxs.first > idxs.second){
        std::swap(idxs.first, idxs.second);
    }

    return idxs;
}




