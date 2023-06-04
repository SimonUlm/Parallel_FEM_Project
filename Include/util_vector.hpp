#ifndef HPC2_UTIL_VECTOR_HPP
#define HPC2_UTIL_VECTOR_HPP

#include <cassert>

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



	/* Vector */
    template<typename T>
    class Vector {
    /*
     * Simple generic vector class which stores elements consectuive in memory for fast
     * memory access
     *
     */
     
    public:
        using ValueType = T;
        using Iterator = SimpleForwardIterator<Vector<T>>;

        /*
    	 * Empty Constructor
         *
    	 */
        explicit Vector() :
                count_(0), data_(nullptr) {
        }
		
		/*
    	 * Constructor to allocate memory for count elements
    	 *
    	 * count: number of elements to be stored
         *
    	 */
        explicit Vector(long count) :
                count_(count), data_(nullptr) {
            if (count != 0)
                data_ = new T[count];
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
            assert(index < count_ || index == 0);
            return data_[index];
        }

        T &operator()(long index) {
            assert(index < count_ || index == 0);
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
        Iterator begin() {
            return Iterator(data_);
        }

        Iterator end() {
            return Iterator(data_ + count_);
        }

    protected:
        long count_;
        T *data_;
    };



    /* BlasVector */
    class BlasVector : public Vector<double> {
        /*
         * Vector for Blas operations.
         *
         * This class implements various Blas Level 1 operations
         *
         */
    public:
        explicit BlasVector(long count) :
                Vector::Vector(count) {}

        // x' * x
        double Dot(BlasVector &y);

        // x <- y
        void Copy(BlasVector &y);

        // x <- alpha * x
        void Scal(double alpha);

        // y <- alpha * x + y
        void Axpy(double alpha, BlasVector &x);

        // max(x)
        double Amax();
    };
}

#endif //HPC2_UTIL_VECTOR_HPP
