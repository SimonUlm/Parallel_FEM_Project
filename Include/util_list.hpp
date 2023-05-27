#ifndef HPC2_MESH_LIST_HPP
#define HPC2_MESH_LIST_HPP

#include "hpc.hpp"


namespace Util {

    template<typename List>
    class ListIterator {
    public:
        using ValueType = typename List::ValueType;

        ListIterator(ValueType *ptr) :
                ptr_(ptr) {}

        ListIterator &operator++() {
            ++ptr_;
            return *this;
        }

        ValueType &operator*() {
            return *ptr_;
        }

        bool operator!=(const ListIterator &other) const {
            return ptr_ != other.ptr_;
        }

    private:
        ValueType *ptr_;
    };

    template<typename T>
    class List {
    public:
        using ValueType = T;
        using Iterator = ListIterator<List<T>>;

        long count;

        explicit List() :
                count(0), data_(nullptr) {
        }

        explicit List(long count) :
                count(count), data_(nullptr) {
            if (count != 0)
                data_ = new T[count];
        }

        ~List() {
            delete[] data_;
        }

        List(List &&other) noexcept:
                data_(other.data_), count(other.count) {
            other.data_ = nullptr;
            other.count = 0;
        }

        List(const List &) = delete;

        List &operator=(List &&other) noexcept {
            delete[] data_;
            data_ = other.data_;
            count = other.count;
            other.data_ = nullptr;
            other.count = 0;
            return *this;
        }

        List &operator=(const List &) = delete;

        T &operator()(long index) const {
            assert(index < count || index == 0);
            return data_[index];
        }

        T &operator()(long index) {
            assert(index < count || index == 0);
            return data_[index];
        }

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
