#include "hpc.hpp"

int main(int argc, char **argv) {
    Util::GeMatrix matrix(4, 5, Util::StorageOrder::ROWMAJOR);
    
    matrix.Init();
    matrix.Print();
    
    Util::SedMatrix sed(5, 25);
    
    sed.Init();
    sed.Print();
    
    Util::GeMatrix check(sed, false);
    
    check.Print();
    
}
