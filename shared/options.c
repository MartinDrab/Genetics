
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "err.h"
#include "utils.h"
#include "options.h"


/************************************************************************/
/*                     GLOBAL VARAIBLES                                 */
/************************************************************************/

static PPROGRAM_OPTION _optionTable = NULL;
static size_t _optionTableSize = 0;
static uint32_t _primes[OPTION_MAX_RECOMMENDED_LENGTH];

/************************************************************************/
/*                    HELPER MACROS                                     */
/************************************************************************/

#define OPTION_GET_VALUE_FUNCTION_INTERNAL(aValueType, aValueTypeName, aOptionType) \
	static ERR_VALUE _option_get_##aValueTypeName(PPROGRAM_OPTION OptionRecord, aValueType *Result) \
	{																		\
		ERR_VALUE ret = ERR_INTERNAL_ERROR;									\
																			\
		assert(flag_on(OptionRecord->Flags, PROGRAM_OPTION_FLAG_IN_USE));	\
		if (OptionRecord->Type == aOptionType) {							\
			*Result = OptionRecord->Value.aValueTypeName;					\
			ret = ERR_SUCCESS;												\
								} else ret = ERR_TYPE_MISMATCH;										\
																			\
		return ret;															\
	}																		\


#define OPTION_GET_VALUE_FUNCTION(aValueType, aValueTypeName, aOptionType) \
	ERR_VALUE option_get_##aValueTypeName(char *OptionName, aValueType *Result) \
	{																			\
		ERR_VALUE ret = ERR_INTERNAL_ERROR;										\
		PPROGRAM_OPTION record = NULL;											\
																				\
		record = _get_option_record(OptionName);								\
		if (record != NULL)														\
			ret = _option_get_##aValueTypeName(record, Result);					\
																																																																else ret = ERR_NOT_FOUND;										\
																				\
		return ret;																\
	}																			\

#define OPTION_ADD_FUNCTION(aValueType, aValueTypeName, aOptionType)					\
	ERR_VALUE option_add_##aValueTypeName(char *OptionName, aValueType DefaultValue)		\
	{																					\
		ERR_VALUE ret = ERR_INTERNAL_ERROR;												\
																						\
		ret = _option_add(OptionName, aOptionType, &DefaultValue, sizeof(aValueType));	\
																						\
		return ret;																		\
	}																					\

#define OPTION_SET_FUNCTION(aValueType, aValueTypeName, aOptionType)					\
	ERR_VALUE option_set_##aValueTypeName(char *OptionName, aValueType DefaultValue)	\
	{																					\
		ERR_VALUE ret = ERR_INTERNAL_ERROR;												\
																						\
		ret = _option_set(OptionName, aOptionType, &DefaultValue);						\
																						\
		return ret;																		\
	}																					\

/************************************************************************/
/*                        HELPER FUNCTIONS                              */
/************************************************************************/


static void _init_primes(void)
{
	size_t primeCount = 0;
	uint32_t currCandidate = 3;

	while (primeCount < OPTION_MAX_RECOMMENDED_LENGTH) {
		if (utils_is_prime(currCandidate)) {
			_primes[primeCount] = currCandidate;
			++primeCount;
		}

		currCandidate += 2;
	}

	return;
}

static size_t _hash_first(char *OptionName)
{
	size_t ret = 0;
	size_t primeIndex = 0;

	while (*OptionName != '\0') {
		ret += (_primes[primeIndex] * *OptionName);
		++OptionName;
		++primeIndex;
		if (primeIndex == OPTION_MAX_RECOMMENDED_LENGTH)
			primeIndex = 0;
	}

	return ret;
}

static size_t _hash_next(const size_t Hash, const size_t Attempt)
{
	assert(Attempt > 0);
	return (Hash + 2 * Attempt + 1);
}

static PPROGRAM_OPTION _get_option_record(char *OptionName)
{
	size_t firstHash = 0;
	PPROGRAM_OPTION ret = NULL;

	firstHash = _hash_first(OptionName) % _optionTableSize;
	ret = _optionTable + firstHash;
	if (flag_on(ret->Flags, PROGRAM_OPTION_FLAG_IN_USE)) {
		if (!strings_equal(ret->Name, OptionName)) {
			size_t attempt = 1;
			size_t hash = firstHash;

			do {
				hash = _hash_next(hash, attempt) % _optionTableSize;
				if (hash == firstHash) {
					ret = NULL;
					break;
				}

				ret = _optionTable + hash;
				if (!flag_on(ret->Flags, PROGRAM_OPTION_FLAG_IN_USE)) {
					ret = NULL;
					break;
				}
				
				++attempt;
			} while (!strings_equal(ret->Name, OptionName));
		}
	} else ret = NULL;

	return ret;
}

static PPROGRAM_OPTION _get_empty_record(char *OptionName)
{
	PPROGRAM_OPTION ret = NULL;
	size_t firstHash = 0;

	firstHash = _hash_first(OptionName) % _optionTableSize;
	ret = _optionTable + firstHash;
	if (flag_on(ret->Flags, PROGRAM_OPTION_FLAG_IN_USE)) {
		size_t attempt = 1;
		size_t hash = firstHash;

		do {
			hash = _hash_next(hash, attempt) % _optionTableSize;
			if (hash == firstHash) {
				ret = NULL;
				break;
			}

			ret = _optionTable + hash;
			++attempt;
		} while (flag_on(ret->Flags, PROGRAM_OPTION_FLAG_IN_USE));
	}

	return ret;
}

static ERR_VALUE _option_add(char *OptionName, EOptionType OptionType, void *DefaultValue, size_t ValueSize)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PPROGRAM_OPTION record = NULL;

	record = _get_empty_record(OptionName);
	if (record != NULL) {
		flag_set(record->Flags, PROGRAM_OPTION_FLAG_IN_USE);
		ret = utils_copy_string(OptionName, &record->Name);
		if (ret == ERR_SUCCESS) {
			record->Type = OptionType;
			if (OptionType != otString) {
				memcpy(&record->Value, DefaultValue, ValueSize);
				record->ValueSize = ValueSize;
				ret = ERR_SUCCESS;
			} else {
				ret = utils_copy_string(*(char **)DefaultValue, &record->Value.String);
				if (ret == ERR_SUCCESS)
					flag_set(record->Flags, PROGRAM_OPTION_FLAG_ALLOCATED);
			}

			if (ret == ERR_SUCCESS)
				flag_set(record->Flags, PROGRAM_OPTION_FLAG_DEFAULT_VALUE);

			if (ret != ERR_SUCCESS)
				utils_free_string(record->Name);
		}

		if (ret != ERR_SUCCESS)
			flag_clear(record->Flags, PROGRAM_OPTION_FLAG_IN_USE);
	} else ret = ERR_TABLE_FULL;

	return ret;
}

static ERR_VALUE _option_set(char *OptionName, EOptionType OptionType, void *Value)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PPROGRAM_OPTION record = NULL;

	record = _get_option_record(OptionName);
	if (record != NULL) {
		if (record->Type == OptionType) {
			if (OptionType != otString) {
				memcpy(&record->Value, Value, record->ValueSize);
				ret = ERR_SUCCESS;
			} else {
				char *tmp = NULL;

				ret = utils_copy_string(*(char **)Value, &tmp);
				if (ret == ERR_SUCCESS) {
					if (flag_on(record->Flags, PROGRAM_OPTION_FLAG_ALLOCATED))
						utils_free_string(record->Value.String);

					record->Value.String = tmp;
					flag_set(record->Flags, PROGRAM_OPTION_FLAG_ALLOCATED);
				}
			}

			if (ret == ERR_SUCCESS)
				flag_clear(record->Flags, PROGRAM_OPTION_FLAG_DEFAULT_VALUE);
		} else ret = ERR_TYPE_MISMATCH;
	} else ret = ERR_NOT_FOUND;

	return ret;
}

static void _option_destroy(PPROGRAM_OPTION Record)
{
	assert(flag_on(Record->Flags, PROGRAM_OPTION_FLAG_IN_USE));
	free(Record->Name);
	Record->Name = NULL;
	if (flag_on(Record->Flags, PROGRAM_OPTION_FLAG_ALLOCATED)) {
		assert(Record->Type == otString);
		free(Record->Value.String);
	}

	Record->Flags = 0;

	return;
}

static ERR_VALUE _set_record_value_str(PPROGRAM_OPTION Record, char *StrValue)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	switch (Record->Type) {
		case otUInt8:
		case otUInt16:
		case otUInt32: {
			char *endptr = StrValue;
			unsigned long result = 0;

			result = strtoul(StrValue, &endptr, 0);
			ret = (result != 0 || endptr != StrValue) ? ERR_SUCCESS : ERR_NOT_AN_UNSIGNED_INTEGRAL_NUMBER;
			if (ret == ERR_SUCCESS) {
				memcpy(&Record->Value, &result, Record->ValueSize);
				flag_clear(Record->Flags, PROGRAM_OPTION_FLAG_DEFAULT_VALUE);
			}
		} break;
		case otInt8:
		case otInt16:
		case otInt32: {
			char *endptr = StrValue;
			long result = 0;

			result = strtol(StrValue, &endptr, 0);
			ret = (result != 0 || endptr != StrValue) ? ERR_SUCCESS : ERR_NOT_AN_INTEGRAL_NUMBER;
			if (ret == ERR_SUCCESS) {
				memcpy(&Record->Value, &result, Record->ValueSize);
				flag_clear(Record->Flags, PROGRAM_OPTION_FLAG_DEFAULT_VALUE);
			}
		} break;
		case otFloat: {
			char *endptr = StrValue;
			float result = 0.0;

			result = strtof(StrValue, &endptr);
			ret = (result != 0 || endptr != StrValue) ? ERR_SUCCESS : ERR_NOT_A_FLOAT_NUMBER;
			if (ret == ERR_SUCCESS) {
				memcpy(&Record->Value, &result, Record->ValueSize);
				flag_clear(Record->Flags, PROGRAM_OPTION_FLAG_DEFAULT_VALUE);
			}
		} break;
		case otDouble: {
			char *endptr = StrValue;
			double result = 0.0;

			result = strtod(StrValue, &endptr);
			ret = (result != 0 || endptr != StrValue) ? ERR_SUCCESS : ERR_NOT_A_DOUBLE_NUMBER;
			if (ret == ERR_SUCCESS) {
				memcpy(&Record->Value, &result, Record->ValueSize);
				flag_clear(Record->Flags, PROGRAM_OPTION_FLAG_DEFAULT_VALUE);
			}
		} break;
		case otString: {
			if (*StrValue != '\0') {
				char *tmp = NULL;

				ret = utils_copy_string(StrValue, &tmp);
				if (ret == ERR_SUCCESS) {
					if (flag_on(Record->Flags, PROGRAM_OPTION_FLAG_ALLOCATED))
						utils_free_string(Record->Value.String);

					Record->Value.String = tmp;
					flag_set(Record->Flags, PROGRAM_OPTION_FLAG_ALLOCATED);
					flag_clear(Record->Flags, PROGRAM_OPTION_FLAG_DEFAULT_VALUE);
				}
			} else ret = ERR_EMPTY_STRING_NOT_ALLOWED;
		} break;
		case otChar: {
			if (strlen(StrValue) == 1) {
				Record->Value.Char = *StrValue;
				flag_clear(Record->Flags, PROGRAM_OPTION_FLAG_DEFAULT_VALUE);
				ret = ERR_SUCCESS;
			} else ret = ERR_STRING_TOO_LONG;
		} break;
		case otBoolean: {
			if (StrValue == NULL || strings_equal(StrValue, "true") || strings_equal(StrValue, "false")) {
				Record->Value.Boolean = (StrValue == NULL || strings_equal(StrValue, "true"));
				flag_clear(Record->Flags, PROGRAM_OPTION_FLAG_DEFAULT_VALUE);
				ret = ERR_SUCCESS;
			} else ret = ERR_INVALID_BOOLEAN_VALUE;
		} break;
		case otUInt64: {
		} break;
		case otInt64: {
		} break;
	}

	return ret;
}

OPTION_GET_VALUE_FUNCTION_INTERNAL(int8_t, Int8, otInt8)
OPTION_GET_VALUE_FUNCTION_INTERNAL(uint8_t, UInt8, otUInt8)
OPTION_GET_VALUE_FUNCTION_INTERNAL(int16_t, Int16, otInt16)
OPTION_GET_VALUE_FUNCTION_INTERNAL(uint16_t, UInt16, otUInt16)
OPTION_GET_VALUE_FUNCTION_INTERNAL(int32_t, Int32, otInt32)
OPTION_GET_VALUE_FUNCTION_INTERNAL(uint32_t, UInt32, otUInt32)
OPTION_GET_VALUE_FUNCTION_INTERNAL(int64_t, Int64, otInt64)
OPTION_GET_VALUE_FUNCTION_INTERNAL(uint64_t, UInt64, otUInt64)
OPTION_GET_VALUE_FUNCTION_INTERNAL(float, Float, otFloat)
OPTION_GET_VALUE_FUNCTION_INTERNAL(double, Double, otDouble)
OPTION_GET_VALUE_FUNCTION_INTERNAL(char, Char, otChar)
OPTION_GET_VALUE_FUNCTION_INTERNAL(char *, String, otString)
OPTION_GET_VALUE_FUNCTION_INTERNAL(boolean, Boolean, otBoolean)


/************************************************************************/
/*                       PUBLIC FUNCTIONS                               */
/************************************************************************/

OPTION_GET_VALUE_FUNCTION(int8_t, Int8, otInt8)
OPTION_GET_VALUE_FUNCTION(uint8_t, UInt8, otUInt8)
OPTION_GET_VALUE_FUNCTION(int16_t, Int16, otInt16)
OPTION_GET_VALUE_FUNCTION(uint16_t, UInt16, otUInt16)
OPTION_GET_VALUE_FUNCTION(int32_t, Int32, otInt32)
OPTION_GET_VALUE_FUNCTION(uint32_t, UInt32, otUInt32)
OPTION_GET_VALUE_FUNCTION(int64_t, Int64, otInt64)
OPTION_GET_VALUE_FUNCTION(uint64_t, UInt64, otUInt64)
OPTION_GET_VALUE_FUNCTION(float, Float, otFloat)
OPTION_GET_VALUE_FUNCTION(double, Double, otDouble)
OPTION_GET_VALUE_FUNCTION(char, Char, otChar)
OPTION_GET_VALUE_FUNCTION(char *, String, otString)
OPTION_GET_VALUE_FUNCTION(boolean, Boolean, otBoolean)

OPTION_ADD_FUNCTION(int8_t, Int8, otInt8)
OPTION_ADD_FUNCTION(uint8_t, UInt8, otUInt8)
OPTION_ADD_FUNCTION(int16_t, Int16, otInt16)
OPTION_ADD_FUNCTION(uint16_t, UInt16, otUInt16)
OPTION_ADD_FUNCTION(int32_t, Int32, otInt32)
OPTION_ADD_FUNCTION(uint32_t, UInt32, otUInt32)
OPTION_ADD_FUNCTION(int64_t, Int64, otInt64)
OPTION_ADD_FUNCTION(uint64_t, UInt64, otUInt64)
OPTION_ADD_FUNCTION(float, Float, otFloat)
OPTION_ADD_FUNCTION(double, Double, otDouble)
OPTION_ADD_FUNCTION(char, Char, otChar)
OPTION_ADD_FUNCTION(char *, String, otString)
OPTION_ADD_FUNCTION(boolean, Boolean, otBoolean)

OPTION_SET_FUNCTION(int8_t, Int8, otInt8)
OPTION_SET_FUNCTION(uint8_t, UInt8, otUInt8)
OPTION_SET_FUNCTION(int16_t, Int16, otInt16)
OPTION_SET_FUNCTION(uint16_t, UInt16, otUInt16)
OPTION_SET_FUNCTION(int32_t, Int32, otInt32)
OPTION_SET_FUNCTION(uint32_t, UInt32, otUInt32)
OPTION_SET_FUNCTION(int64_t, Int64, otInt64)
OPTION_SET_FUNCTION(uint64_t, UInt64, otUInt64)
OPTION_SET_FUNCTION(float, Float, otFloat)
OPTION_SET_FUNCTION(double, Double, otDouble)
OPTION_SET_FUNCTION(char, Char, otChar)
OPTION_SET_FUNCTION(char *, String, otString)
OPTION_SET_FUNCTION(boolean, Boolean, otBoolean)


ERR_VALUE options_parse_command_line(int argc, char **argv)
{
	int i = 0;
	char *argName = NULL;
	size_t argnameLen = 0;
	char *argValue = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	while (ret == ERR_SUCCESS && i < argc) {
		argName = argv[i];
		argnameLen = strlen(argName);
		argValue = NULL;
		if (argnameLen > 2 && memcmp(argName, "--", 2 * sizeof(char)) == 0) {
			PPROGRAM_OPTION record = NULL;
			
			argName += 2;
			argnameLen -= 2;
			record = _get_empty_record(argName);
			if (record != NULL) {
				if (record->Type != otBoolean) {
					if (i < argc - 1) {
						++i;
						argValue = argv[i];
					} else ret = ERR_OPTION_VALUE_NOT_FOUND;
				}

				if (ret == ERR_SUCCESS)
					ret = _set_record_value_str(record, argValue);
			} else ret = ERR_UNKNOWN_OPTION;
		}

		++i;
	}

	return ret;
}


/************************************************************************/
/*                      INITIALIZATION AND FINALIZATION                 */
/************************************************************************/

ERR_VALUE options_module_init(const size_t MaxOptions)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (utils_is_prime(MaxOptions)) {
		_optionTableSize = MaxOptions;
		_optionTable = (PPROGRAM_OPTION)calloc(_optionTableSize, sizeof(PROGRAM_OPTION));
		if (_optionTable != NULL) {
			_init_primes();
			ret = ERR_SUCCESS;
		} else ret = ERR_OUT_OF_MEMORY;
	} else ret = ERR_NOT_A_PRIME;

	return ret;
}

void options_module_finit(void)
{
	PPROGRAM_OPTION record = _optionTable;

	for (size_t i = 0; i < _optionTableSize; ++i) {
		if (flag_on(record->Flags, PROGRAM_OPTION_FLAG_IN_USE))
			_option_destroy(record);

		++record;
	}

	return;
}
