#ifndef HPC2_MESH_OBJECTS_HPP
#define HPC2_MESH_OBJECTS_HPP

#include <algorithm>
#include <cstdio>

#include "hpc.hpp"


namespace Mesh{
    struct Node {
        double x; // x coordinate
        double y; // y coordinate

        // Add node to node coordinate-wise
        Node operator+(const Node other) const {
            return Node{x + other.x, y + other.y};
        }

        // Subtract node from node coordinate-wise
        Node operator-(const Node other) const {
            return Node{x - other.x, y - other.y};
        }

        // Multiply node coords by scalar
        Node operator*(const double mult) const {
            return Node{mult * x, mult * y};
        }

        // Dot product of two nodes interpreted as vectors
        double dot(Node n) {
        	return x * n.x + y * n.y;
        }
        
        // Print node data
        void Print();
    };

    struct Element {
        long n1, n2, n3;  // nodes
        long m1, m2, m3;  // edges
        long t;           // affiliation

        // Utility functions for mesh refinement
        long get_n(long n);

        long get_successor_n(long n);
        
        long get_predecessor_n(long n); 

        long get_m(long m); 

        long get_successor_m(long m); 

        long get_predecessor_m(long m); 
        
        // Print element data
        void Print();
    };

    struct Edge {
        long n1; // starting node
        long n2; // end node
    };

    struct BoundaryEdge {
        long n1, n2; // nodes
        long m;      // edge
        long t;      // 0 if Neumann boundary, 1 if Dirichlet boundary
        
        // Print boundary data
        void Print();
    };
    
}

#endif //HPC2_MESH_OBJECtS_HPP
