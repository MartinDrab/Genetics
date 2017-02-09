
#ifdef WIN32
#include <Windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#endif
#include "err.h"
#include "utils.h"
#include "file-utils.h"




ERR_VALUE utils_file_map(const char *FileName, PFUTILS_MAPPED_FILE Handle)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
#ifdef WIN32
	HANDLE mapHandle = NULL;
	HANDLE fileHandle = INVALID_HANDLE_VALUE;
#else
	int fd;
	struct stat st;
#endif

#ifdef WIN32
	fileHandle = CreateFileA(FileName, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (fileHandle != INVALID_HANDLE_VALUE) {
		LARGE_INTEGER fileSize;

		if (GetFileSizeEx(fileHandle, &fileSize)) {
			Handle->Size = fileSize.QuadPart;
			mapHandle = CreateFileMapping(fileHandle, NULL, PAGE_READONLY, 0, 0, NULL);
			if (mapHandle != NULL) {
				Handle->Address = MapViewOfFile(mapHandle, FILE_MAP_READ, 0, 0, 0);
				if (Handle->Address != NULL)
					ret = ERR_SUCCESS;

				CloseHandle(mapHandle);
			}
			else ret = ERR_IO_ERROR;
		}
		else ret = ERR_IO_ERROR;

		CloseHandle(fileHandle);
	} else ret = ERR_NOT_FOUND;
#else
	stat(FileName, &st);
	fd = open(FileName, O_RDONLY);
	if (fd != -1) {
		Handle->Size = st.st_size;
		Handle->Address = mmap(0, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
		if (Handle->Address != MAP_FAILED)
			ret = ERR_SUCCESS;
	
		close(fd);
	} else ret = ERR_IO_ERROR;
#endif

	return ret;
}


void utils_file_unmap(PFUTILS_MAPPED_FILE Handle)
{

#ifdef WIN32
	UnmapViewOfFile(Handle->Address);
#else
	munmap(Handle->Address, (size_t)Handle->Size);
#endif

	return;
}




ERR_VALUE utils_fopen(const char *FileName, const uint32_t Mode, FILE **Stream)
{
	FILE *tmpStream = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char *strModes[] = {
		"",
		"rb",
		"wb",
		"wb+"
		"ab",
		"rb",
		"wb",
		"wb+",
		"",
		"r",
		"w",
		"w+"
		"a",
		"r",
		"w",
		"w+"
	};

#pragma warning(disable : 4996)
	tmpStream = fopen(FileName, strModes[Mode]);
	if (tmpStream != NULL) {
		*Stream = tmpStream;
		ret = ERR_SUCCESS;
	}
	else ret = ERR_IO_ERROR;

	return ret;
}


ERR_VALUE utils_fread(void *Buffer, const size_t Size, const size_t Count, FILE *Stream)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	if (fread(Buffer, Size, Count, Stream) != Count)
		ret = ERR_IO_ERROR;

	return ret;
}


ERR_VALUE utils_fwrite(const void *Buffer, const size_t Size, const size_t Count, FILE *Stream)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	if (fwrite(Buffer, Size, Count, Stream) != Count)
		ret = ERR_IO_ERROR;

	return ret;
}


ERR_VALUE utils_fclose(FILE *Stream)
{
	return (fclose(Stream) == 0) ? ERR_SUCCESS : ERR_IO_ERROR;
}


ERR_VALUE utils_file_read(const char *FileName, char **Data, size_t *DataLength)
{
	FILE *f = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_fopen(FileName, FOPEN_MODE_READ, &f);
	if (f != NULL) {
#ifdef _MSC_VER
		ret = _fseeki64(f, 0, SEEK_END);
#else
		ret = fseeko(f, 0, SEEK_END);
#endif
		if (ret == ERR_SUCCESS) {
			off_t fileSize = 0;
#ifdef _MSC_VER
			fileSize = _ftelli64(f);
#else
			fileSize = ftello(f);
#endif
			if (fileSize != -1L) {
#ifdef _MSC_VER
				ret = _fseeki64(f, 0, SEEK_SET);
#else
				ret = fseeko(f, 0, SEEK_SET);
#endif
				if (ret == ERR_SUCCESS) {
					char *tmpData = NULL;
					size_t tmpSize = (size_t)fileSize;

					ret = utils_malloc(tmpSize + sizeof(char), &tmpData);
					if (ret == ERR_SUCCESS) {
						ret = utils_fread(tmpData, 1, tmpSize, f);
						if (ret == ERR_SUCCESS) {
							tmpData[tmpSize] = 0;
							*Data = tmpData;
							*DataLength = tmpSize;
						}

						if (ret != ERR_SUCCESS)
							utils_free(tmpData);
					}
				} else ret = ERR_IO_ERROR;
			} else ret = ERR_IO_ERROR;
		} else ret = ERR_IO_ERROR;

		utils_fclose(f);
	}

	return ret;
}
