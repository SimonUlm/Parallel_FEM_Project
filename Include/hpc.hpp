#ifndef HPC2_HPC_HPP
#define HPC2_HPC_HPP

#include "mesh_objects.hpp"
#include "util_list.hpp"
#include "vector_converter.hpp"
#include "mesh.hpp"
#include "skeleton.hpp"
#include "parallel_dot_product.hpp"

#ifdef _MPI
void MpiPrintSerial(Mesh::LocalMesh &mesh, Skeleton::Skeleton &skeleton, MPI_Comm comm, int rank, int nof_processes);
void MpiPrintVectorSerial(std::vector<long> &vector_1, std::vector<long> &vector_2, MPI_Comm comm, int rank, int nof_processes);
#endif

#endif //HPC2_HPC_HPP
