/* bambi_utils.h

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

#ifndef BAMBI_UTILS_H
#define BAMBI_UTILS_H

void store_msg(char **str, const char *fmt, ...);
void display(const char *fmt, ...);
void die(const char *fmt, ...);

#define smalloc(s) _s_malloc((s), __FILE__, __LINE__, __func__)
#define srealloc(p, s) _s_realloc((p), (s), __FILE__, __LINE__, __func__)
void * _s_malloc(size_t size, const char *file, unsigned int line, const char *func);
void * _s_realloc(void *ptr, size_t size, const char *file, unsigned int line, const char *func);


#endif
