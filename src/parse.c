/*  parse.c -- simple parsing for arguments

    Copyright (C) 2016 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "parse.h"

/*
 * Parse a comma delimited set of arguments into an array
 */
void parse_tags(va_t *tags, char *arg)
{
    char *argstr = strdup(arg);
    char *s = strtok(argstr,",");
    while (s) {
        va_push(tags,strdup(s));
        s = strtok(NULL,",");
    }
    free(argstr);
}

void parse_int(ia_t *ia, char *arg)
{
    char *argstr = strdup(arg);
    char *s = strtok(argstr,",");
    while (s) {
        ia_push(ia,atoi(s));
        s = strtok(NULL,",");
    }
    free(argstr);
}

