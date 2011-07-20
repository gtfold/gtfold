#include "utils.h"
#include <string.h>
#include <string>

char baseToDigit(const char* base) {
	if (!strcmp(base, "A")) {
		return '1';
	}
	if (!strcmp(base, "C")) {
		return '2';
	}
	if (!strcmp(base, "G")) {
		return '3';
	}
	if (!strcmp(base, "U")) {
		return '4';
	}
	if (!strcmp(base, "N")) {
		return '5';
	}
	return (char) '0';
}

/**
 * Encodes the given base to a number.
 *
 * Handles upper- and lower-case letters.
 *
 * @param base The character to encode
 * @return The numerical encoding of the base
 */
unsigned char encode(char base) {
	switch(base) {
		case 'A':
		case 'a':
			return BASE_A;
		case 'C':
		case 'c':
			return BASE_C;
		case 'G':
		case 'g':
			return BASE_G;
		case 'U':
		case 'u':
		case 'T':
		case 't':
			return BASE_U;
		case 'N':
		case 'n':
		case 'R':
		case 'r':
		case 'Y':
		case 'y':
		case 'M':
		case 'm':
		case 'K':
		case 'k':
		case 'S':
		case 's':
		case 'W':
		case 'w':
		case 'B':
		case 'b':
		case 'D':
		case 'd':
		case 'H':
		case 'h':
		case 'V':
		case 'v':
			return BASE_A;
		default:
			return 'X';
	}
}

/**
 * Returns true if the given base is an A, C, G, U, or T
 * @param base 'A', 'C', 'G', 'U', etc.
 */
int isWatsonCrickBase(char base) {
	switch(base) {
		case 'A':
		case 'a':
		case 'C':
		case 'c':
		case 'G':
		case 'g':
		case 'U':
		case 'u':
		case 'T':
		case 't':
			return true;
		default:
			return false;
	}
}
