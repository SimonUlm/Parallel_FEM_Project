#ifndef HPC2_MESH_HPP
#define HPC2_MESH_HPP

#include <cassert>

using namespace std;

namespace Mesh{

    struct Node {
        double x;
        double y;
    };

    struct Element {
        long n1, n2, n3;
        long m1, m2, m3;
        long t;
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

    template<typename T>
    class List {
    public:
        long count;

        explicit List(long count) :
                count(count), data(new T[count]) {
        }
        ~List() {
            delete [] data;
        }
        List(const List &) = delete;

        T & operator()(long index) const {
            assert(index < count);
            return data[index];
        }

        T & operator()(long index) {
            assert(index < count);
            return data[index];
        }

    private:
        T *data;
    };

    class RectangularMesh {
    private:
        long m;
        long n;

    public:
        List<Node> nodes;
        List<Element> elements;
        List<Edge> edges;
        List<BoundaryEdge> boundary;
        List<long> fixed_nodes;

        RectangularMesh(long m, long n) :
                m(m), n(n),
                nodes((m + 1) * (n + 1)), elements(2 * m * n), edges(3 * m * n + m + n),
                boundary(2 * (m + n)), fixed_nodes(m + n + 1) {
        }

        RectangularMesh(long nnodes, long nelem, long nedges, long nbdry, long nfixed) :
                m(1), n(1),
                nodes(nnodes), elements(nelem), edges(nedges),
                boundary(nbdry), fixed_nodes(nfixed) {
        }

        // Defined in create.cpp
        void create();

        // TODO
        //void refine();
    };
}

#endif //HPC2_MESH_HPP
