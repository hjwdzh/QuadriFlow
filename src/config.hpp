#ifndef CONFIG_H_
#define CONFIG_H_

// Move settings to cmake to make git happy :)

// #define WITH_SCALE
// #define WITH_CUDA
// #define WITH_TBB
// #define WITH_OPENGL

#define DATA_PATH "./data"
#define INVALID -1
#define GRAIN_SIZE 1024

#define LOG_OUTPUT

#include <chrono>

// simulation of Windows GetTickCount()
unsigned long long inline GetCurrentTime64() {
    using namespace std::chrono;
    return duration_cast<milliseconds>(steady_clock::now().time_since_epoch()).count();
}
#endif
