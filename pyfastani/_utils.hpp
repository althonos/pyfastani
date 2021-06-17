#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <stdint.h>
#include <chrono>
#include <limits>
#include <iostream>
#include <vector>
#include <zlib.h>

#include "common/kseq.h"
#include "map/include/base_types.hpp"
#include "map/include/winSketch.hpp"

// not needed anywhere except in `cgi::correctRefGenomeIds`
// so we can just patch them
extern int omp_get_thread_num(void);
extern int omp_get_num_threads(void);

typedef kseq_t* kseq_ptr_t;

static const char COMPLEMENT_LOOKUP[128] = {
    '\x00', '\x01', '\x02', '\x03', '\x04', '\x05', '\x06', '\x07',
    '\x08', '\t',   '\n',   '\x0', '\x0c', '\r',   '\x0e', '\x0f',
    '\x10', '\x11', '\x12', '\x13', '\x14', '\x15', '\x16', '\x17',
    '\x18', '\x19', '\x1a', '\x1', '\x1c', '\x1d', '\x1e', '\x1f',
    ' ',    '!',    '"',    '#',    '$',    '%',    '&',    '\'',
    '(',    ')',    '*',    '+',    ',',    '-',    '.',    '/',
    '0',    '1',    '2',    '3',    '4',    '5',    '6',    '7',
    '8',    '9',    ':',    ';',    '<',    '=',    '>',    '?',
    '@',    'T',    'V',    'G',    'H',    'E',    'F',    'C',
    'D',    'I',    'J',    'M',    'L',    'K',    'N',    'O',
    'P',    'Q',    'Y',    'S',    'A',    'U',    'B',    'W',
    'X',    'R',    'Z',    '[',    '\\',   ']',    '^',    '_',
    '`',    't',    'v',    'g',    'h',    'e',    'f',    'c',
    'd',    'i',    'j',    'm',    'l',    'k',    'n',    'o',
    'p',    'q',    'y',    's',    'a',    'u',    'b',    'w',
    'x',    'r',    'z',    '{',    '|',    '}',    '~',    '\x7f'
};

inline char complement(char base) {
    return COMPLEMENT_LOOKUP[(size_t) (base & 0x7F)];
}

#endif
