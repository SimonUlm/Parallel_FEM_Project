#ifndef HPC2_VECTOR_CONVERTER_HPP
#define HPC2_VECTOR_CONVERTER_HPP

#ifdef _MPI
#include <mpi.h>
#endif

namespace Conversion {

    class VectorConverter {
    public:
        VectorConverter() :
                n_global_nodes_(0),
                local_nodes_priority_(),
                local_to_global_(nullptr) {}

        VectorConverter(Util::List<long> global_nodes_priority, Util::List<long> &local_to_global) :
                n_global_nodes_(global_nodes_priority.count),
                local_nodes_priority_(local_to_global.count),
                local_to_global_(&local_to_global) {
            for (long i = 0; i < local_nodes_priority_.count; ++i)
                local_nodes_priority_(i) = global_nodes_priority(local_to_global(i));
        }

        VectorConverter(VectorConverter &&) = delete;
        VectorConverter(const VectorConverter &) = delete;

        VectorConverter & operator=(VectorConverter &&other) = default;
        VectorConverter & operator=(const VectorConverter &) = delete;

#ifdef _MPI
        void AccumulatedToDistributed(std::vector<long> &vector);
        void DistributedToAccumulated(std::vector<long> &vector, MPI_Comm comm);
#endif

    private:
        long n_global_nodes_;
        Util::List<long> local_nodes_priority_;
        Util::List<long> *local_to_global_;
    };
}

#endif //HPC2_VECTOR_CONVERTER_HPP
