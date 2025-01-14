#ifndef HPC2_UTIL_VECTOR_HPP
#define HPC2_UTIL_VECTOR_HPP

#ifndef NDEBUG

#include <cassert>

#endif

#include <iostream>

#include "hpc.hpp"


namespace Util {

    template<typename List>
    class SimpleIterator {
    public:
        using ValueType = typename List::ValueType;
        using PointerType = ValueType *;
        using ReferenceType = ValueType &;

        explicit SimpleIterator(ValueType *ptr) :
                ptr_(ptr) {}

        SimpleIterator &operator++() {
            ++ptr_;
            return *this;
        }

        const SimpleIterator operator++(int) {
            SimpleIterator iterator = *this;
            ++(*this);
            return *this;
        }

        SimpleIterator &operator--() {
            --ptr_;
            return *this;
        }

        const SimpleIterator operator--(int) {
            SimpleIterator iterator = *this;
            --(*this);
            return *this;
        }

        ReferenceType operator[](int index) {
            return ptr_[index];
        }

        PointerType operator->() {
            return ptr_;
        }

        ReferenceType operator*() {
            return *ptr_;
        }

        bool operator!=(const SimpleIterator &other) const {
            return ptr_ != other.ptr_;
        }

        bool operator==(const SimpleIterator &other) const {
            return !(*this == other);
        }

    private:
        ValueType *ptr_;
    };


    /*
     * Vector
     *
     * Simple generic vector class which stores elements consectuive in memory for fast
     * memory access
     *
     */
    template<typename T>
    class Vector {
    public:
        using ValueType = T;
        using Iterator = SimpleIterator<Vector<T>>;

        /*
    	 * Default Constructor
         *
    	 */
        explicit Vector() :
                count_(0), data_(nullptr) {
        }

        /*
         * Constructor to allocate memory for count elements
         *
         * count_: number of elements to be stored
         *
         */
        explicit Vector(long count) :
                count_(count), data_(nullptr) {
            if (count != 0)
                data_ = new T[count]{};
        }

        // Destructor
        ~Vector() {
            delete[] data_;
        }

        Vector(Vector &&other) noexcept:
                data_(other.data_), count_(other.count_) {
            other.data_ = nullptr;
            other.count_ = 0;
        }

        Vector(const Vector &) = delete;

        // Assignment operations
        Vector &operator=(Vector &&other) noexcept {
            delete[] data_;
            data_ = other.data_;
            count_ = other.count_;
            other.data_ = nullptr;
            other.count_ = 0;
            return *this;
        }

        Vector &operator=(const Vector &) = delete;

        // Access operations and getters
        const T &operator()(long index) const {
#ifndef NDEBUG
            assert(index < count_ || index == 0);
#endif
            return data_[index];
        }

        T &operator()(long index) {
#ifndef NDEBUG
            assert(index < count_ || index == 0);
#endif
            return data_[index];
        }

        const long count() const { return count_; }

        T *const data() const { return data_; }

        // Initializers
        void Init() {
            for (long i = 0; i < count_; ++i)
                data_[i] = i + 1;
        }

        void Init(T value) {
            for (long i = 0; i < count_; ++i)
                data_[i] = value;
        }

        // Iterators
        const Iterator begin() const {
            return Iterator(data_);
        }

        Iterator begin() {
            return Iterator(data_);
        }

        const Iterator end() const {
            return Iterator(data_ + count_);
        }

        Iterator end() {
            return Iterator(data_ + count_);
        }

        // Other
        void Copy(Vector &other) {
#ifndef NDEBUG
            assert(count() == other.count());
#endif
            for (long i = 0; i < count_; ++i)
                data_[i] = other(i);
        }

        void Print() const {
            std::cout << std::endl << "Vector:" << std::endl;
            for (long i = 0; i < count_; ++i)
                std::cout << data_[i] << std::endl;
            std::cout << std::endl;
        }

    protected:
        long count_;
        T *data_;
    };


    /*
     * BlasVector
     *
     * Vector for Blas operations
     *
     * This class implements various Blas Level 1 operations
     *
     */
    class BlasVector : public Vector<double> {
    public:
        explicit BlasVector() :
                Vector::Vector() {}

        explicit BlasVector(long count) :
                Vector::Vector(count) {}

        BlasVector(BlasVector &&other) noexcept:
                Vector::Vector(std::move(other)) {}

        BlasVector(const BlasVector &) = delete;

        // Assignment operators
        BlasVector &operator=(BlasVector &&other) noexcept {
            Vector::operator=(std::move(other));
            return *this;
        }

        BlasVector &operator=(const BlasVector &) = delete;

        // x' * x
        double Dot(BlasVector &y);

        // x <- alpha * x
        void Scal(double alpha);

        // y <- alpha * x + y
        void Axpy(double alpha, BlasVector &x);

        // max(x)
        double Amax();
    };
}

#endif //HPC2_UTIL_VECTOR_HPP
