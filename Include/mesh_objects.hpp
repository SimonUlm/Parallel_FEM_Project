#ifndef HPC2_MESH_OBJECTS_HPP
#define HPC2_MESH_OBJECTS_HPP

#include <algorithm>
#include <cassert>


namespace Mesh{
    struct Node {
        double x;
        double y;

        Node operator+(const Node other) const {
            return Node{x + other.x, y + other.y};
        }

        Node operator*(const double mult) const {
            return Node{mult * x, mult * y};
        }
    };

    struct Element {
        long n1, n2, n3;
        long m1, m2, m3;
        long t;

        long get_n(long n) const {
            switch (n) {
                case 1:
                    return n1;
                case 2:
                    return n2;
                case 3:
                    return n3;
                default:
                    return n1;
            }
        }

        long get_successor_n(long n) const {
            n = (n == 3) ? 0 : n;
            return get_n(n+1);
        }

        long get_predecessor_n(long n) const {
            n = (n == 1) ? 4 : n;
            return get_n(n-1);
        }

        long get_m(long m) const {
            switch (m) {
                case 1:
                    return m1;
                case 2:
                    return m2;
                case 3:
                    return m3;
                default:
                    return m1;
            }
        }

        long get_successor_m(long m) const {
            m = (m == 3) ? 0 : m;
            return get_m(m+1);
        }

        long get_predecessor_m(long m) const {
            m = (m == 1) ? 4 : m;
            return get_m(m-1);
        }
    };

    struct Edge {
        long n1;
        long n2;
    };

    struct BoundaryEdge {
        long n1, n2;
        long m;
        long t;
    };
    
}

#endif //HPC2_MESH_OBJECtS_HPP
