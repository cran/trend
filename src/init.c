#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "trend.h"
/* 
    Copyright (C) 2017-2018 Thorsten Pohlert

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* define argument types */
static R_NativePrimitiveArgType mcbr_f[] = {
  REALSXP, INTSXP, INTSXP, REALSXP};

/* define Fortran entry points, their names, nr of arguments*/
static const R_FortranMethodDef FortEntries[]  = {
  {"mcbr", (DL_FUNC) &F77_SUB(mcbr), 4, mcbr_f},
  {"mcbu", (DL_FUNC) &F77_SUB(mcbu), 4, mcbr_f},
  {"mcsnht", (DL_FUNC) &F77_SUB(mcsnht), 4, mcbr_f},
  {NULL, 0}
};

// register
void attribute_visible R_init_trend(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
