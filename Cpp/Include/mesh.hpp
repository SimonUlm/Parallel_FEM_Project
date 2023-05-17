#ifndef HPC2_MESH_HPP
#define HPC2_MESH_HPP

#include <algorithm>
#include <cassert>
#include <mpi.h>



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
            switch (n) {
                case 1:
                    return n2;
                case 2:
                    return n3;
                case 3:
                    return n1;
                default:
                    return n1;
            }
        }

        long get_predecessor_n(long n) const {
            switch (n) {
                case 1:
                    return n3;
                case 2:
                    return n1;
                case 3:
                    return n2;
                default:
                    return n1;
            }
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
            switch (m) {
                case 1:
                    return m2;
                case 2:
                    return m3;
                case 3:
                    return m1;
                default:
                    return m1;
            }
        }

        long get_predecessor_m(long m) const {
            switch (m) {
                case 1:
                    return m3;
                case 2:
                    return m1;
                case 3:
                    return m2;
                default:
                    return m1;
            }
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

    template<typename T>
    class List {
    private:
        T *data;

    public:
        long count;

        explicit List() :
                count(0), data(nullptr) {
        }
        explicit List(long count) :
                count(count), data(nullptr) {
            if (count != 0)
                data = new T[count];
        }
        ~List() {
            delete [] data;
        }
        List(List &&) = delete;
        List(const List &) = delete;

        List & operator=(List &&other) noexcept {
            delete [] data;
            data = other.data;
            count = other.count;
            other.data = nullptr;
            other.count = 0;
            return *this;
        }
        List & operator=(const List &) = delete;

        T & operator()(long index) const {
            assert(index < count || index == 0);
            return data[index];
        }
        T & operator()(long index) {
            assert(index < count || index == 0);
            return data[index];
        }
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

        RectangularMesh() :
                m(0), n(0) {
        }
        RectangularMesh(long m, long n) :
                m(m), n(n),
                nodes((m + 1) * (n + 1)), elements(2 * m * n), edges(3 * m * n + m + n),
                boundary(2 * (m + n)), fixed_nodes(m + n + 1) {
        }
        RectangularMesh(long m, long n, long nnodes, long nelem, long nbdry) :
                m(m), n(n),
                nodes(nnodes), elements(nelem), edges(0),
                boundary(nbdry), fixed_nodes(0) {
        }
        RectangularMesh(RectangularMesh &&) = delete;
        RectangularMesh(const RectangularMesh &) = delete;

        RectangularMesh & operator=(RectangularMesh &&other) = default;
        RectangularMesh & operator=(const RectangularMesh &) = delete;

        // Defined in Create.cpp
        void Create(Node bottom_left_node = Node{0, 0}, Node top_right_node = Node{1, 1});

        // Defined in Refine.cpp
        void Refine();

        // Defined in Print.cpp
        void Print();

        // Defined in Scatter.cpp
        void Scatter(RectangularMesh &local_mesh, MPI_Comm comm, int rank, int nof_local_elem);
    };
}



// Other declarations
void MpiPrintSerial(Mesh::RectangularMesh &mesh, MPI_Comm comm, int rank, int nof_processes);

#endif //HPC2_MESH_HPP
