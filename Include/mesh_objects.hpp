#ifndef HPC2_MESH_OBJECTS_HPP
#define HPC2_MESH_OBJECTS_HPP

#include <algorithm>
#include <cassert>
#include <cstdio>

#include "hpc.hpp"


namespace Mesh{
    struct Node {
        double x; // x-Coordinate of Node
        double y; // y-Coordinate of Node

        // Add coordinates component-wise
        Node operator+(const Node other) const {
            return Node{x + other.x, y + other.y};
        }

        // Subtract coordinates component-wise
        Node operator-(const Node other) const {
            return Node{x - other.x, y - other.y};
        }

        // Scale coordinates
        Node operator*(const double mult) const {
            return Node{mult * x, mult * y};
        }

        // Calculate dot product of the coordinates of two Nodes
        double dot(Node n) {
        	return x * n.x + y * n.y;
        }
        
        void Print(); // Print x and y coordinates of the Node
    };

    struct Element {
        long n1, n2, n3; // Node Numbers of the Element
        long m1, m2, m3; // Edge Numbers of the Element
        long t;          // Affiliation of the Element

        long get_n(long n);             // Get Node n

        long get_successor_n(long n);   // Get Node n+1, where the successor of Node 3 is 1
        
        long get_predecessor_n(long n); // GetNode n-1, where the predecessor of Node 1 is 3

        long get_m(long m);             // Get Edge m

        long get_successor_m(long m);   // Get Edge n+1, where the successor of Edge 3 is 1

        long get_predecessor_m(long m); // Get Edge n-1, where the predecessor of Edge 1 is 3
        
        void Print();       	        // Print all Nodes all Edges and the Affiliation belonging to the Element
    };

    struct Edge {
        long n1;    // Start Node of the Edge
        long n2;    // End Node of the Edge
    };

    struct BoundaryEdge {
        long n1, n2;    // Start and end Node of the Edge
        long m;         // Edge number
        long t;         // Type of the Edge (0: Neumann, 1: Dirichlet)
        
        void Print();   // Print Nodes, Edge numbers and Type of the Boundary Edge
    };
    
}

#endif //HPC2_MESH_OBJECtS_HPP
