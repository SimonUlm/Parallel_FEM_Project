#include <cstdio>

#include "hpc.hpp"

namespace Mesh {
    // Get Node n
    const long Element::get_n(long n) const {
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

    // Get Node n+1, where the successor of Node 3 is 1
    const long Element::get_successor_n(long n) const {
        n = (n == 3) ? 0 : n;
        return get_n(n + 1);
    }

    // Get Node n-1, where the predecessor of Node 1 is 3
    const long Element::get_predecessor_n(long n) const {
        n = (n == 1) ? 4 : n;
        return get_n(n - 1);
    }

    // Get Edge m
    const long Element::get_m(long m) const {
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

    // Get Edge m+1, where the successor of Edge 3 is 1
    const long Element::get_successor_m(long m) const {
        m = (m == 3) ? 0 : m;
        return get_m(m + 1);
    }

    const long Element::get_predecessor_m(long m) const {
        // Get Edge m-1, where the predecessor of Edge 1 is 3
        m = (m == 1) ? 4 : m;
        return get_m(m - 1);
    }
}
