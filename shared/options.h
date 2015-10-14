
#ifndef __PROGRAM_OPTIONS_H__
#define __PROGRAM_OPTIONS_H__

#include <stdint.h>
#include "err.h"
#include "utils.h"


#define OPTION_MAX_RECOMMENDED_LENGTH			64

typedef enum _EOptionType {
	otUnknown,
	otInt8,
	otUInt8,
	otInt16,
	otUInt16,
	otInt32,
	otUInt32,
	otInt64,
	otUInt64,
	otFloat,
	otDouble,
	otChar,
	otString,
	otBoolean,
	otMaximumType,
} EOptionType, *PEOptionType;

#define PROGRAM_OPTION_FLAG_IN_USE					0x1
#define PROGRAM_OPTION_FLAG_ALLOCATED				0x2
#define PROGRAM_OPTION_FLAG_DEFAULT_VALUE			0x4
#define PROGRAM_OPTION_FLAG_DESCRIPTION_ALLOCATED	0x8

typedef struct _PROGRAM_PTION {
	uint8_t Flags;
	char Shortcut;
	char *Name;
	char *Description;
	EOptionType Type;
	size_t ValueSize;
	union _PROGRAM_OPTION_VALUE_DATA {
		int8_t Int8;
		uint8_t UInt8;
		int16_t Int16;
		uint16_t UInt16;
		int32_t Int32;
		uint32_t UInt32;
		int64_t Int64;
		uint64_t UInt64;
		float Float;
		double Double;
		char Char;
		char *String;
		boolean Boolean;
	} Value;
	unsigned int Order;
} PROGRAM_OPTION, *PPROGRAM_OPTION;


#define OPTION_GET_VALUE_FUNCTION_HEADER(aValueType, aValueTypeName) \
	ERR_VALUE option_get_##aValueTypeName(const char *OptionName, aValueType *Result)

OPTION_GET_VALUE_FUNCTION_HEADER(int8_t, Int8);
OPTION_GET_VALUE_FUNCTION_HEADER(uint8_t, UInt8);
OPTION_GET_VALUE_FUNCTION_HEADER(int16_t, Int16);
OPTION_GET_VALUE_FUNCTION_HEADER(uint16_t, UInt16);
OPTION_GET_VALUE_FUNCTION_HEADER(int32_t, Int32);
OPTION_GET_VALUE_FUNCTION_HEADER(uint32_t, UInt32);
OPTION_GET_VALUE_FUNCTION_HEADER(int64_t, Int64);
OPTION_GET_VALUE_FUNCTION_HEADER(uint64_t, UInt64);
OPTION_GET_VALUE_FUNCTION_HEADER(float, Float);
OPTION_GET_VALUE_FUNCTION_HEADER(double, Double);
OPTION_GET_VALUE_FUNCTION_HEADER(char, Char);
OPTION_GET_VALUE_FUNCTION_HEADER(char *, String);
OPTION_GET_VALUE_FUNCTION_HEADER(boolean, Boolean);

#define OPTION_ADD_FUNCTION_HEADER(aValueType, aValueTypeName)	\
	ERR_VALUE option_add_##aValueTypeName(const char *OptionName, const aValueType DefaultValue)	

OPTION_ADD_FUNCTION_HEADER(int8_t, Int8);
OPTION_ADD_FUNCTION_HEADER(uint8_t, UInt8);
OPTION_ADD_FUNCTION_HEADER(int16_t, Int16);
OPTION_ADD_FUNCTION_HEADER(uint16_t, UInt16);
OPTION_ADD_FUNCTION_HEADER(int32_t, Int32);
OPTION_ADD_FUNCTION_HEADER(uint32_t, UInt32);
OPTION_ADD_FUNCTION_HEADER(int64_t, Int64);
OPTION_ADD_FUNCTION_HEADER(uint64_t, UInt64);
OPTION_ADD_FUNCTION_HEADER(float, Float);
OPTION_ADD_FUNCTION_HEADER(double, Double);
OPTION_ADD_FUNCTION_HEADER(char, Char);
OPTION_ADD_FUNCTION_HEADER(char *, String);
OPTION_ADD_FUNCTION_HEADER(boolean, Boolean);

#define OPTION_SET_FUNCTION_HEADER(aValueType, aValueTypeName)	\
	ERR_VALUE option_set_##aValueTypeName(const char *OptionName, const aValueType DefaultValue)	

OPTION_SET_FUNCTION_HEADER(int8_t, Int8);
OPTION_SET_FUNCTION_HEADER(uint8_t, UInt8);
OPTION_SET_FUNCTION_HEADER(int16_t, Int16);
OPTION_SET_FUNCTION_HEADER(uint16_t, UInt16);
OPTION_SET_FUNCTION_HEADER(int32_t, Int32);
OPTION_SET_FUNCTION_HEADER(uint32_t, UInt32);
OPTION_SET_FUNCTION_HEADER(int64_t, Int64);
OPTION_SET_FUNCTION_HEADER(uint64_t, UInt64);
OPTION_SET_FUNCTION_HEADER(float, Float);
OPTION_SET_FUNCTION_HEADER(double, Double);
OPTION_SET_FUNCTION_HEADER(char, Char);
OPTION_SET_FUNCTION_HEADER(char *, String);
OPTION_SET_FUNCTION_HEADER(boolean, Boolean);

ERR_VALUE option_set_description_const(const char *Name, const char *Description);
ERR_VALUE opttion_set_description(const char *Name, const char *Description);
ERR_VALUE option_set_shortcut(const char *Name, const char Shortcut);


ERR_VALUE options_parse_command_line(int argc, char **argv);
void options_print(void);
void options_print_help(void);



ERR_VALUE options_module_init(const size_t MaxOptions);
void options_module_finit(void);


#endif 
