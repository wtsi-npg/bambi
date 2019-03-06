/*  bambi_utils.c -- bambi utility functions

    Copyright (C) 2016 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#include "bambi_utils.h"

void store_msg(char **str, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    free(*str);
    *str = smalloc(1024);
    vsnprintf(*str, 1023, fmt, ap);
}

void display(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
}

void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap,fmt);
    fflush(stdout);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

void * _s_malloc(size_t size, const char *file, unsigned int line, const char *func)
{
    void *m = malloc(size);
    if (!m) die("Couldn't allocate %zd bytes in %s at %s line %u: %s\n", size, func, file, line, strerror(errno));
    return m;
}

void * _s_realloc(void *ptr, size_t size, const char *file, unsigned int line, const char *func)
{
    void *m = realloc(ptr, size);
    if (!m) die("Couldn't reallocate %zd bytes in %s at %s line %u: %s\n", size, func, file, line, strerror(errno));
    return m;
}
