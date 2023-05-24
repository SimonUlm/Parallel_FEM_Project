#ifndef HPC2_HPC_HPP
#define HPC2_HPC_HPP

#include "mesh_objects.hpp"
#include "util_list.hpp"
#include "mesh.hpp"
#include "skeleton.hpp"

#ifdef _MPI
void MpiPrintSerial(Mesh::LocalMesh &mesh, Skeleton::Skeleton &skeleton, MPI_Comm comm, int rank, int nof_processes);
#endif

#endif //HPC2_HPC_HPP
