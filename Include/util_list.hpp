#ifndef HPC2_MESH_LIST_HPP
#define HPC2_MESH_LIST_HPP

#include "hpc.hpp"


namespace Util {
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
}

#endif //HPC2_MESH_LIST_HPP
