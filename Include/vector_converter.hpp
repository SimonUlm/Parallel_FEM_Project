#ifndef HPC2_VECTOR_CONVERTER_HPP
#define HPC2_VECTOR_CONVERTER_HPP

#ifdef _MPI
#include <mpi.h>
#endif

namespace Util {

    /*!  \class VectorConverter vector_converter_.hpp "Include/vector_converter_.hpp"
     *   \brief Converts different vector types into each other
     *
     *	 The VectorConverter class is used for converting accumulated vectors into distributed vectors and vice versa
     */
    class VectorConverter {
    public:
        VectorConverter() :
                n_global_nodes_(0),
                local_nodes_priority_(),
                local_to_global_(nullptr) {}

        VectorConverter(Vector<long> &global_nodes_priority, Vector<long> &local_to_global) :
                n_global_nodes_(global_nodes_priority.count()),
                local_nodes_priority_(local_to_global.count()),
                local_to_global_(&local_to_global) {
            for (long i = 0; i < local_nodes_priority_.count(); ++i)
                local_nodes_priority_(i) = global_nodes_priority(local_to_global(i));
        }

        VectorConverter(VectorConverter &&) = delete;
        VectorConverter(const VectorConverter &) = delete;

        VectorConverter & operator=(VectorConverter &&other) = default;
        VectorConverter & operator=(const VectorConverter &) = delete;

#ifdef _MPI
        /*!
    	 *   Converts accumulated vector into distributed vector
    	 */
        void AccumulatedToDistributed(Vector<double> &vector) const;

        /*!
    	 *   Converts distributed vector into accumulated vector
    	 */
        void DistributedToAccumulated(Vector<double> &vector, MPI_Comm comm) const;
#endif

    private:
        long n_global_nodes_; /*!< number of global nodes */
        Vector<long> local_nodes_priority_; /*!< counts how many processes share each node  */
        Vector<long> *local_to_global_; /*!< reference to vector that maps local to global nodes */
    };
}

#endif //HPC2_VECTOR_CONVERTER_HPP
