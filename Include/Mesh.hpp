#ifndef HPC2_MESH_HPP
#define HPC2_MESH_HPP

#include <cassert>

namespace Meshes{

    class Mesh {
        // Declare data types
    private:
        struct Node {
            double x;
            double y;
        };
        
        struct Element {
            int n1;
            int n2;
            int n3;
            int m1;
            int m2;
            int m3;
            int t;
        };

        struct Edge {
            int n1;
            int n2;
        };

        struct Boundary {
            int n1;
            int n2;
            int m;
            int t;
        };

        template<typename T>
        class List {
        public:
            int count;

            explicit List(int count) :
                    count(count), data(new T[count]) {
            }
            ~List() {
                delete [] data;
            }
            List(const List &) = delete;

            T & operator()(int i) {
                assert(i < count);
                return data[i];
            }
        private:
            T *data;
        };

        // Defines actual class
    public:
        List<Node> nodes;
        List<Element> elements;
        List<Edge> edges;
        List<Boundary> boundaries;
        List<int> fixed;

        Mesh(int nnodes, int nelem, int nedges, int nbdry, int nfixed) :
            nodes(nnodes), elements(nelem), edges(nedges),
            boundaries(nbdry), fixed(nfixed) {
        }
    };
}

#endif //HPC2_MESH_HPP
