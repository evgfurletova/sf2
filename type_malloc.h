#ifndef TYPE_MALLOC
#define TYPE_MALLOC

#include <stdlib.h>
#include <string.h>

template<typename T>
inline T * Malloc(size_t NumberOfElements)
{
    void * ptr = malloc(NumberOfElements * sizeof(T));
    memset(ptr, 0x00, NumberOfElements * sizeof(T));
    T * result = reinterpret_cast<T*>(ptr);
    return result;
}

/*
 * Usage:
 *         int * blablabla = Malloc<int>(100);
 *         double * aaa = Malloc<double>(50);
 * 
 *         free(blablabla);
 *         free(aaa);
 * 
 */
#endif
