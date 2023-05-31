#ifndef HPC2_MESH_LIST_HPP
#define HPC2_MESH_LIST_HPP

#include "hpc.hpp"


namespace Util {

    template<typename List>
    class SimpleForwardIterator {
    public:
        using ValueType = typename List::ValueType;

        explicit SimpleForwardIterator(ValueType *ptr) :
                ptr_(ptr) {}

        SimpleForwardIterator &operator++() {
            ++ptr_;
            return *this;
        }

        ValueType &operator*() {
            return *ptr_;
        }

        bool operator!=(const SimpleForwardIterator &other) const {
            return ptr_ != other.ptr_;
        }

    private:
        ValueType *ptr_;
    };

    template<typename T>
    class Vector {
    public:
        using ValueType = T;
        using Iterator = SimpleForwardIterator<Vector<T>>;

        long count;

        // Constructors
        explicit Vector() :
                count(0), data_(nullptr) {
        }

        explicit Vector(long count) :
                count(count), data_(nullptr) {
            if (count != 0)
                data_ = new T[count];
        }

        ~Vector() {
            delete[] data_;
        }

        Vector(Vector &&other) noexcept:
                data_(other.data_), count(other.count) {
            other.data_ = nullptr;
            other.count = 0;
        }

        Vector(const Vector &) = delete;

        // Assignment operations
        Vector &operator=(Vector &&other) noexcept {
            delete[] data_;
            data_ = other.data_;
            count = other.count;
            other.data_ = nullptr;
            other.count = 0;
            return *this;
        }

        Vector &operator=(const Vector &) = delete;

        // Access operations and getter
        const T & operator()(long index) const {
            assert(index < count || index == 0);
            return data_[index];
        }

        T & operator()(long index) {
            assert(index < count || index == 0);
            return data_[index];
        }

        // Iterators
        Iterator begin() {
            return Iterator(data_);
        }

        Iterator end() {
            return Iterator(data_ + count);
        }

    private:
        T *data_;
    };
}

#endif //HPC2_MESH_LIST_HPP
