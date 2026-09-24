// memory.c
#if defined(__APPLE__)
#include <mach/mach.h>
#elif defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#elif defined(__linux__)
#include <sys/time.h>
#include <sys/resource.h>
#endif

long long get_peak_memory_bytes() {
#if defined(__APPLE__)
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) != KERN_SUCCESS) {
        return -1;
    }
    return (long long)info.resident_size_max;

#elif defined(_WIN32)
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc))) {
        return (long long)pmc.PeakWorkingSetSize;
    }
    return -1;

#elif defined(__linux__)
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        return (long long)usage.ru_maxrss * 1024;
    }
    return -1;

#else
    return -1;
#endif
}
