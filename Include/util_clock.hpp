#ifndef CLOCK_HPP
#define CLOCK_HPP

#include <time.h>
#include <sys/time.h>

namespace Util {
    double get_wall_time();
    double get_cpu_time();
}

#endif // CLOCK_HPP 