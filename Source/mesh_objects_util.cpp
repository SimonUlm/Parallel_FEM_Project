#include <cstdio>

#include "hpc.hpp"

namespace Mesh {
    long Element::get_n(long n) {
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

    long Element::get_successor_n(long n) {
    	n = (n == 3) ? 0 : n;
        return get_n(n+1);
    }

    long Element::get_predecessor_n(long n) {
        n = (n == 1) ? 4 : n;
        return get_n(n-1);
    }

    long Element::get_m(long m) {
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

    long Element::get_successor_m(long m) {
    	m = (m == 3) ? 0 : m;
        return get_m(m+1);
    }

    long Element::get_predecessor_m(long m) {
    	m = (m == 1) ? 4 : m;
        return get_m(m-1);
    }
}
