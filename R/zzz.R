## Copyright (C) 2022-2023        Ching-Chuan Chen
##
## This file is part of oneMKLUtil.
##
## oneMKLUtil is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## oneMKLUtil is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with oneMKLUtil. If not, see <http://www.gnu.org/licenses/>.

#' @importFrom rlang env_unlock env_binding_unlock env_binding_lock env_lock

old_base_solve <- NULL

.onAttach <- function(libname, pkgname) {
  # replace solve due to the broken issue
  old_base_solve <- base::solve
  packageStartupMessage("Replace base::solve with fMatSolve due to the bug to import MKL routines.\n")
  rlang::env_unlock(env = asNamespace('base'))
  rlang::env_binding_unlock(env = asNamespace('base'))
  assign('solve', fMatSolve, envir = asNamespace('base'))
  rlang::env_binding_lock(env = asNamespace('base'))
  rlang::env_lock(asNamespace('base'))
}

.onDetach <- function(libpath) {
  packageStartupMessage("Recover base::solve.\n")
  rlang::env_unlock(env = asNamespace('base'))
  rlang::env_binding_unlock(env = asNamespace('base'))
  assign('solve', old_base_solve, envir = asNamespace('base'))
  rlang::env_binding_lock(env = asNamespace('base'))
  rlang::env_lock(asNamespace('base'))
}
