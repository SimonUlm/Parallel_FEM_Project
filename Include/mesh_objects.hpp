#ifndef HPC2_MESH_OBJECTS_HPP
#define HPC2_MESH_OBJECTS_HPP

#include <algorithm>
#include <cassert>
#include <cstdio>

#include "hpc.hpp"


namespace Mesh{
    struct Node {
        double x;
        double y;

        Node operator+(const Node other) const {
            return Node{x + other.x, y + other.y};
        }
        
        Node operator-(const Node other) const {
            return Node{x - other.x, y - other.y};
        }

        Node operator*(const double mult) const {
            return Node{mult * x, mult * y};
        }
        
        double dot(Node n) {
        	return x * n.x + y * n.y;
        }
        
        void Print();       
    };

    struct Element {
        long n1, n2, n3;
        long m1, m2, m3;
        long t;

        long get_n(long n);

        long get_successor_n(long n);
        
        long get_predecessor_n(long n); 

        long get_m(long m); 

        long get_successor_m(long m); 

        long get_predecessor_m(long m); 
        
        void Print();
    };

    struct Edge {
        long n1;
        long n2;
    };

    struct BoundaryEdge {
        long n1, n2;
        long m;
        long t; //  0 if Neumann boundary, 1 if Dirichlet boundary
        
        void Print();
    };
    
}

#endif //HPC2_MESH_OBJECtS_HPP
