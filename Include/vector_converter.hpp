#ifndef HPC2_VECTOR_CONVERTER_HPP
#define HPC2_VECTOR_CONVERTER_HPP

#ifdef _MPI
#include <mpi.h>
#endif

namespace Skeleton {

    // Declare Skeleton class
    class Skeleton;

    /*!  \class VectorConverter vector_converter_.hpp "Include/vector_converter_.hpp"
     *   \brief Converts different vector types into each other
     *
     *	 The VectorConverter class is used for converting accumulated vectors into distributed vectors and vice versa
     */
    class VectorConverter {
    public:
        VectorConverter() :
                n_global_nodes_(0),
                n_global_crosspoints_(0),
                local_nodes_priority_(),
                local_to_global_(nullptr) {}

        VectorConverter(long m, long n,
                        Util::Vector<long> &global_nodes_priority, Util::Vector<long> &local_to_global) :
                n_global_nodes_(global_nodes_priority.count()),
                n_global_crosspoints_((m + 1) * (n + 1)),
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
        /*
    	 *   Converts accumulated vector into distributed vector (in-place)
    	 */
        void AccumulatedToDistributed(Util::Vector<double> &local_vector) const;

        /*
    	 *   Converts distributed vector into accumulated vector (in-place)
    	 */
        void DistributedToAccumulated(Util::Vector<double> &local_vector,
                                      const Skeleton &local_skel) const;

        /*
    	 *   Takes distributed vector as input and returns accumulated vector
    	 */
        void DistributedToAccumulated(Util::Vector<double> &local_vector_send,
                                      Util::Vector<double> &local_vector_recv,
                                      const Skeleton &local_skel) const;

        /*
         * Gather accumulated local vector and store into global vector
         */
        void GatherAccumulatedVector(Util::Vector<double> &local_vector_send,
                                     Util::Vector<double> &global_vector_recv,
                                     const Skeleton &local_skel) const;

        /*
         * Gather accumulated local vector and store into global vector
         */
        void GatherDistributedVector(Util::Vector<double> &local_vector_send,
                                     Util::Vector<double> &global_vector_recv,
                                     const Skeleton &local_skel) const;
#endif

    private:
        long n_global_nodes_;
        long n_global_crosspoints_; // number of global nodes
        Util::Vector<long> local_nodes_priority_; // counts how many processes share each node
        Util::Vector<long> *local_to_global_; // reference to vector that maps local to global nodes
    };
}

#endif //HPC2_VECTOR_CONVERTER_HPP
