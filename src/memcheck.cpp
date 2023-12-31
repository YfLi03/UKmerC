/* Taken from https://hpcf.umbc.edu/general-productivity/checking-memory-usage/ */
#include "memcheck.hpp"
#include <upcxx/upcxx.hpp>

int get_memory_usage_kb(long* vmrss_kb, long* vmsize_kb)
{
    /* Get the the current process' status file from the proc filesystem */
    FILE* procfile = fopen("/proc/self/status", "r");

    long to_read = 8192;
    char buffer[to_read];
    int read = fread(buffer, sizeof(char), to_read, procfile);
    fclose(procfile);

    short found_vmrss = 0;
    short found_vmsize = 0;
    char* search_result;

    /* Look through proc status contents line by line */
    char delims[] = "\n";
    char* line = strtok(buffer, delims);

    while (line != NULL && (found_vmrss == 0 || found_vmsize == 0) )
    {
        search_result = strstr(line, "VmRSS:");
        if (search_result != NULL)
        {
            sscanf(line, "%*s %ld", vmrss_kb);
            found_vmrss = 1;
        }

        search_result = strstr(line, "VmSize:");
        if (search_result != NULL)
        {
            sscanf(line, "%*s %ld", vmsize_kb);
            found_vmsize = 1;
        }

        line = strtok(NULL, delims);
    }

    return (found_vmrss == 1 && found_vmsize == 1) ? 0 : 1;
}

int get_cluster_memory_usage_kb(size_t& vmrss, size_t& vmsize, int root, int np)
{
    long vmrss_kb;
    long vmsize_kb;
    int ret_code = get_memory_usage_kb(&vmrss_kb, &vmsize_kb);

    if (ret_code != 0)
    {
        printf("Could not gather memory usage!\n");
        return ret_code;
    }

    upcxx::reduce_one(vmrss_kb, upcxx::op_fast_add, root).wait();
    upcxx::reduce_one(vmsize_kb, upcxx::op_fast_add, root).wait();

    vmrss = vmrss_kb;
    vmsize = vmsize_kb;

    return 0;
}


void get_mem(int nprocs, int myrank, size_t& vmrss, size_t& vmsize){
    size_t vmrss_kb;
    size_t vmsize_kb;
    get_cluster_memory_usage_kb(vmrss_kb, vmsize_kb, 0, nprocs);


    if (myrank == 0) {
        vmrss = vmrss_kb / (1024 * 1024);
        vmsize = vmsize_kb / (1024 * 1024);
    }
}