
#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"
#include "read-info.h"



void read_info_init(PREAD_INFO Info)
{
	dym_array_init_READ_INFO_ENTRY(&Info->Array, 140);

	return;
}


ERR_VALUE read_info_assign(PREAD_INFO Info, const GEN_ARRAY_READ_INFO_ENTRY *Array)
{
	dym_array_clear_READ_INFO_ENTRY(&Info->Array);
	return dym_array_push_back_array_READ_INFO_ENTRY(&Info->Array, Array);
}


void read_info_finit(PREAD_INFO Info)
{
	dym_array_finit_READ_INFO_ENTRY(&Info->Array);

	return;
}


ERR_VALUE read_info_copy(PREAD_INFO Dest, const READ_INFO *Source)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const size_t count = read_info_get_count(Source);
	const READ_INFO_ENTRY *entry = Source->Array.Data;

	for (size_t i = 0; i < count; ++i) {
		ret = read_info_add(Dest, entry->ReadIndex, entry->ReadPosition);
		if (ret != ERR_SUCCESS)
			break;

		++entry;
	}

	return ret;
}


ERR_VALUE read_info_add(PREAD_INFO Info, const size_t ReadIndex, const size_t ReadPosition)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	READ_INFO_ENTRY entry;

	entry.ReadIndex = ReadIndex;
	entry.ReadPosition = ReadPosition;
	ret = dym_array_push_back_READ_INFO_ENTRY(&Info->Array, entry);

	return ret;
}


ERR_VALUE read_info_intersection(PREAD_INFO Info1, PREAD_INFO Info2, GEN_ARRAY_PTYPE(READ_INFO_ENTRY) Intersection, const boolean AscendingPosition, const size_t ReadDistance)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const READ_INFO_ENTRY *entry1 = Info1->Array.Data;
	const READ_INFO_ENTRY *entry2 = Info2->Array.Data;
	size_t index1 = 0;
	size_t index2 = 0;
	const size_t count1 = read_info_get_count(Info1);
	const size_t count2 = read_info_get_count(Info2);

	ret = ERR_SUCCESS;
	while (ret == ERR_SUCCESS && index1 < count1 && index2 < count2) {
		if (entry1->ReadIndex == entry2->ReadIndex) {
			if (!AscendingPosition || (entry1->ReadPosition <= entry2->ReadPosition && (ReadDistance == 0 || entry1->ReadPosition + ReadDistance == entry2->ReadPosition)))
				ret = dym_array_push_back_READ_INFO_ENTRY(Intersection, *entry1);
			
			++entry1;
			++index1;
			++entry2;
			++index2;
		} else if (entry1->ReadIndex < entry2->ReadIndex) {
			++entry1;
			++index1;
		} else {
			++entry2;
			++index2;
		}
	}

	return ret;
}


ERR_VALUE read_info_diff(const READ_INFO *Info, const GEN_ARRAY_TYPE(READ_INFO_ENTRY) *Subtrahend, GEN_ARRAY_PTYPE(READ_INFO_ENTRY) Difference)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const READ_INFO_ENTRY *entry1 = Info->Array.Data;
	const READ_INFO_ENTRY *entry2 = Subtrahend->Data;
	size_t index1 = 0;
	size_t index2 = 0;
	const size_t count1 = read_info_get_count(Info);
	const size_t count2 = gen_array_size(Subtrahend);

	ret = ERR_SUCCESS;
	while (ret == ERR_SUCCESS && index1 < count1 && index2 < count2) {
		if (entry1->ReadIndex == entry2->ReadIndex) {
			++entry1;
			++index1;
			++entry2;
			++index2;
		} else if (entry1->ReadIndex < entry2->ReadIndex) {
			ret = dym_array_push_back_READ_INFO_ENTRY(Difference, *entry1);
			++entry1;
			++index1;
		} else {
			++entry2;
			++index2;
		}
	}

	while (ret == ERR_SUCCESS && index1 < count1) {
		ret = dym_array_push_back_READ_INFO_ENTRY(Difference, *entry1);
		++entry1;
		++index1;
	}

	return ret;
}
