/* bambi.h

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

#ifndef __BAMBI_H__
#define __BAMBI_H__

#include "config.h"
#include "array.h"
#include "bambi_utils.h"

#define INDEX_SEPARATOR "-"
#define QUAL_SEPARATOR " "

// Machine Type is used by i2b 
typedef enum { MT_UNKNOWN,
               MT_MISEQ,           // MiSeq and HiSeq 2000/2500
               MT_NEXTSEQ,         // MiniSeq and NextSeq 500/550
               MT_HISEQX,          // HiSeq X and HiSeq 3000/4000
               MT_NOVASEQ          // NovaSeq
             } MACHINE_TYPE;

const char *bambi_version(void);

#endif

