
#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"
#include "read-info.h"


double read_info_weight(const READ_INFO *Info, const size_t CurrentReadIndex, const size_t CurrentReadPosition)
{
	double ret = 0;
	size_t readIndex = (size_t)-1;
	const READ_INFO_ENTRY *entry = Info->Array.Data;

	for (size_t i = 0; i < read_info_get_count(Info); ++i) {
		if (readIndex == entry->ReadIndex) {
			++entry;
			continue;
		}
		
		if (entry->ReadIndex == CurrentReadIndex &&
			entry->ReadPosition < CurrentReadPosition) {
			++entry;
			continue;
		}

		if (entry->Quality == 0)
			ret += 0;
		else if (entry->Quality < 20)
			ret += 0.5;
		else if (entry->Quality < 30)
			ret += 0.75;
		else if (entry->Quality < 40)
			ret += 1;
		else if (entry->Quality == 255)
			ret += 20;
		else ret += 1;

		readIndex = entry->ReadIndex;
		++entry;
	}

	return ret;
}

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
	return read_info_add_array(Dest, &Source->Array);
}


ERR_VALUE read_info_add_array(PREAD_INFO Info, const GEN_ARRAY_READ_INFO_ENTRY *Array)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const size_t count = gen_array_size(Array);
	const READ_INFO_ENTRY *entry = Array->Data;

	for (size_t i = 0; i < count; ++i) {
		ret = read_info_add(Info, entry->ReadIndex, entry->ReadPosition, entry->Quality);
		if (ret != ERR_SUCCESS)
			break;

		++entry;
	}

	return ret;
}


ERR_VALUE read_info_add(PREAD_INFO Info, const size_t ReadIndex, const size_t ReadPosition, const uint8_t Quality)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	READ_INFO_ENTRY entry;

	entry.ReadIndex = ReadIndex;
	entry.ReadPosition = ReadPosition;
	entry.Quality = Quality;
	ret = dym_array_push_back_READ_INFO_ENTRY(&Info->Array, entry);

	return ret;
}


void read_info_remove(PREAD_INFO Info, const size_t ReadIndex, const size_t ReadPosition)
{
	READ_INFO_ENTRY entry;
	PREAD_INFO_ENTRY item = Info->Array.Data;

	entry.ReadIndex = ReadIndex;
	entry.ReadPosition = ReadPosition;
	for (size_t i = 0; i < gen_array_size(&Info->Array); ++i) {
		if (item->ReadIndex == entry.ReadIndex &&
			item->ReadPosition == entry.ReadPosition) {
			memmove(item, item + 1, (read_info_get_count(Info) - i - 1)*sizeof(READ_INFO_ENTRY));
			dym_array_pop_back_READ_INFO_ENTRY(&Info->Array);
			break;
		}

		++item;
	}

	return;
}


void read_info_remove_2(PREAD_INFO Info, const size_t ReadIndex)
{
	boolean removed = FALSE;

	do {
		PREAD_INFO_ENTRY entry = Info->Array.Data;

		removed = FALSE;
		for (size_t i = 0; i < read_info_get_count(Info); ++i) {
			removed = (entry->ReadIndex == ReadIndex);
			if (removed) {
				memmove(entry, entry + 1, (read_info_get_count(Info) - i - 1)*sizeof(READ_INFO_ENTRY));
				dym_array_pop_back_READ_INFO_ENTRY(&Info->Array);
				break;
			}

			++entry;
		}
	} while (removed);

	return;
}


ERR_VALUE read_info_intersection(PREAD_INFO Info1, PREAD_INFO Info2, GEN_ARRAY_PTYPE(READ_INFO_ENTRY) Intersection, const size_t ReadDistance)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const READ_INFO_ENTRY *entry1 = Info1->Array.Data;
	const READ_INFO_ENTRY *entry2 = Info2->Array.Data;
	size_t index1 = 0;
	size_t index2 = 0;
	const size_t count1 = read_info_get_count(Info1);
	const size_t count2 = read_info_get_count(Info2);
	size_t readIndex = (size_t)-1;

	ret = ERR_SUCCESS;
	dym_array_clear_READ_INFO_ENTRY(Intersection);
	while (ret == ERR_SUCCESS && index1 < count1 && index2 < count2) {
		if (entry1->ReadIndex == entry2->ReadIndex) {
			if (entry1->ReadPosition + ReadDistance == entry2->ReadPosition) {
				readIndex = entry1->ReadIndex;
				ret = dym_array_push_back_READ_INFO_ENTRY(Intersection, *entry1);
				++entry1;
				++index1;
				++entry2;
				++index2;
			} else {
				if (entry1->ReadPosition + ReadDistance < entry2->ReadPosition) {
					++entry1;
					++index1;
				} else {
					++entry2;
					++index2;
				}
			}
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


void read_info_subtract(PREAD_INFO Info, const GEN_ARRAY_TYPE(READ_INFO_ENTRY) *Subtrahend, const size_t Distance)
{
	const READ_INFO_ENTRY *entry1 = Info->Array.Data;
	const READ_INFO_ENTRY *entry2 = Subtrahend->Data;
	size_t index1 = 0;
	size_t index2 = 0;
	const size_t count1 = read_info_get_count(Info);
	const size_t count2 = gen_array_size(Subtrahend);

	Info->Array.ValidLength = 0;
	while (index1 < count1 && index2 < count2) {
		if (entry1->ReadIndex == entry2->ReadIndex) {
			if (entry1->ReadPosition - Distance == entry2->ReadPosition) {
				++entry1;
				++index1;
				++entry2;
				++index2;
			} else if (entry1->ReadPosition - Distance < entry2->ReadPosition) {
				dym_array_push_back_no_alloc_READ_INFO_ENTRY(&Info->Array, *entry1);
				++entry1;
				++index1;
			} else {
				++entry2;
				++index2;
			}
		} else if (entry1->ReadIndex < entry2->ReadIndex) {
			dym_array_push_back_no_alloc_READ_INFO_ENTRY(&Info->Array, *entry1);
			++entry1;
			++index1;
		} else {
			++entry2;
			++index2;
		}
	}

	while (index1 < count1) {
		dym_array_push_back_no_alloc_READ_INFO_ENTRY(&Info->Array, *entry1);
		++entry1;
		++index1;
	}

	return;
}


static int _read_info_entry_compare(const READ_INFO_ENTRY *E1, const READ_INFO_ENTRY *E2)
{
	int res = ((int)E1->ReadIndex - (int)E2->ReadIndex);

	if (res == 0)
		res = ((int)E1->ReadPosition - (int)E2->ReadPosition);

	return res;
}

ERR_VALUE read_info_union(PREAD_INFO Target, const GEN_ARRAY_READ_INFO_ENTRY *Source)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const size_t count1 = read_info_get_count(Target);
	const size_t count2 = gen_array_size(Source);
	size_t index1 = 0;
	size_t index2 = 0;
	const READ_INFO_ENTRY *entry1 = Target->Array.Data;
	const READ_INFO_ENTRY *entry2 = Source->Data;

	ret = ERR_SUCCESS;
	while (ret == ERR_SUCCESS && index1 < count1 && index2 < count2) {
		if (entry1->ReadIndex == entry2->ReadIndex) {
			if (entry1->ReadPosition < entry2->ReadPosition) {
				++entry1;
				++index1;
			} else if (entry1->ReadPosition == entry2->ReadPosition) {
				++entry2;
				++index2;
			} else {
				ret = dym_array_push_back_READ_INFO_ENTRY(&Target->Array, *entry2);
				++entry2;
				++index2;
			}
		} else if (entry1->ReadIndex < entry2->ReadIndex) {
			++entry1;
			++index1;
		} else {
			ret = dym_array_push_back_READ_INFO_ENTRY(&Target->Array, *entry2);
			++entry2;
			++index2;
		}
	}

	while (ret == ERR_SUCCESS && index2 < count2) {
		ret = dym_array_push_back_READ_INFO_ENTRY(&Target->Array, *entry2);
		++entry2;
		++index2;
	}

	qsort(Target->Array.Data, gen_array_size(&Target->Array), sizeof(READ_INFO_ENTRY), _read_info_entry_compare);

	return ret;
}

ERR_VALUE read_info_merge(PREAD_INFO Dest, const READ_INFO *Info1, const READ_INFO *Info2)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const READ_INFO_ENTRY *entry1 = NULL;
	size_t index1 = 0;
	const READ_INFO_ENTRY *entry2 = NULL;
	size_t index2 = 0;
	const size_t count1 = read_info_get_count(Info1);
	const size_t count2 = read_info_get_count(Info2);

	ret = dym_array_reserve_READ_INFO_ENTRY(&Dest->Array, count1 + count2);
	if (ret == ERR_SUCCESS) {
		size_t lastReadIndex = (size_t)-1;
		
		entry1 = Info1->Array.Data;
		entry2 = Info2->Array.Data;
		while (index1 < count1 && index2 < count2) {
			while (index1 < count1 && entry1->ReadIndex == lastReadIndex) {
				++entry1;
				++index1;
			}
			
			while (index2 < count2 && entry2->ReadIndex == lastReadIndex) {
				++entry2;
				++index2;
			}

			if (index1 < count1 && index2 < count2) {
				if (entry1->ReadIndex <= entry2->ReadIndex) {
					dym_array_push_back_no_alloc_READ_INFO_ENTRY(&Dest->Array, *entry1);
					lastReadIndex = entry1->ReadIndex;
				} else {
					dym_array_push_back_no_alloc_READ_INFO_ENTRY(&Dest->Array, *entry2);
					lastReadIndex = entry2->ReadIndex;
				}
			}
		}

		size_t tmp = lastReadIndex;
		while (index1 < count1) {
			if (entry1->ReadIndex != tmp) {
				dym_array_push_back_no_alloc_READ_INFO_ENTRY(&Dest->Array, *entry1);
				tmp = entry1->ReadIndex;
			}

			++entry1;
			++index1;
		}

		tmp = lastReadIndex;
		while (index2 < count2) {
			if (entry2->ReadIndex != tmp) {
				dym_array_push_back_no_alloc_READ_INFO_ENTRY(&Dest->Array, *entry2);
				tmp = entry2->ReadIndex;
			}

			++entry2;
			++index2;
		}
	}

	return ret;
}
