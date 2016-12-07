
#ifndef __FILE_UTILS_H__
#define __FILE_UTILS_H__


#include "err.h"




#define FOPEN_MODE_READ				1
#define FOPEN_MODE_WRITE			2
#define FOPEN_MODE_APPEND			4

typedef struct _FUILTS_MAPPED_FILE {
	void *Address;
	uint64_t Size;
} FUTILS_MAPPED_FILE, *PFUTILS_MAPPED_FILE;


ERR_VALUE utils_file_map(const char *FileName, PFUTILS_MAPPED_FILE Handle);
void utils_file_unmap(PFUTILS_MAPPED_FILE Handle);

#define utils_in_mapped(aHandle, aPointer)	((uintptr_t)(aHandle.Address) <= (uintptr_t)(aPointer) && (uintptr_t)(aPointer) < (uintptr_t)(aHandle).Address + (uintptr_t)(aHandle).Size)

ERR_VALUE utils_file_read(const char *FileName, char **Data, size_t *DataLength);

ERR_VALUE utils_fopen(const char *FileName, const uint32_t Mode, FILE **Stream);
ERR_VALUE utils_fread(void *Buffer, const size_t Size, const size_t Count, FILE *Stream);
ERR_VALUE utils_fwrite(const void *Buffer, const size_t Size, const size_t Count, FILE *Stream);
ERR_VALUE utils_fclose(FILE *Stream);





#endif
