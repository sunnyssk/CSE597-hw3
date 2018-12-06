#include "meminfo.h"

void GetMemUsage (char * res_container) {
    FILE *fin = fopen("/proc/self/status", "r");
    int result = -1;
    char buffer[128];
    while (fgets(buffer, 128, fin) != NULL) {
        if(strncmp("VmSize:", buffer, 7) == 0) strcpy(res_container, buffer + 7);
    }
    fclose(fin);
}

void MemSizeOutput (char * buffer) {
    GetMemUsage(buffer);
    printf("Current memory usage: %s", buffer);
}