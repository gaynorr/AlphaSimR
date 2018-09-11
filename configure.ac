#                                               -*- Autoconf -*-
# AlphaSimR configure.ac
# (based on RcppArmadillo configure.ac)
# 
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(AlphaSimR, m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]))

## Set R_HOME
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
	AC_MSG_ERROR([Could not determine R_HOME.])
fi

## Use R to set CXX and CXXFLAGS
CXX=$(${R_HOME}/bin/R CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

## We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP

## Check the C++ compiler using the CXX value set
AC_PROG_CXX

use_openmp=""
if test "${GXX}" = yes; then
	AC_MSG_CHECKING([g++ version is sufficient])
	gxx_version=$(${CXX} -v 2>&1 | awk '/^.*g.. version/ {print $3}')
	
	case ${gxx_version} in
		1.*|2.*|3.*|4.0.*|4.1.*|4.2.*|4.3.*|4.4.*|4.5.*|4.6.*|4.7.0|4.7.1)
			AC_MSG_RESULT([no])
			AC_MSG_WARN([Must use g++ version 4.7.2 or greater])
			AC_MSG_ERROR([Use a different compiler])
		;;
		4.7.*|4.8.*|4.9.*|5.0*|5.1*|5.2*|5.3*)
			AC_MSG_RESULT([yes, but no OpenMP])
			use_openmp="no"
		;;
		5.4*|5.5*|5.6*|5.7*|5.8*|5.9*|6.*|7.*|8.*|9.*|10.*)
			AC_MSG_RESULT([yes, with OpenMP])
			use_openmp="yes"
		;;
		*)
			AC_MSG_RESULT([not g++])
			use_openmp="no"
		;;
	esac
fi

## Check for Apple
AC_MSG_CHECKING([for macOS])
RSysinfoName=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')

if test x"${RSysinfoName}" == x"Darwin"; then
	AC_MSG_RESULT([yes])
	AC_MSG_CHECKING([for macOS Apple compiler])
	apple_compiler=$($CXX --version 2>&1 | grep -i -c -e 'apple llvm')
	
	if test x"${apple_compiler}" == x"1"; then
		AC_MSG_RESULT([yes])
		AC_MSG_WARN([OpenMP unavailable])
		use_openmp="no"
	else
		AC_MSG_RESULT([no])
		AC_MSG_CHECKING([for clang compiler])
		clang_compiler=$($CXX --version 2>&1 | grep -i -c -e 'clang ')
		
		if test x"${clang_compiler}" == x"1"; then
			AC_MSG_RESULT([yes])
			AC_MSG_CHECKING([for OpenMP compatible])
			clang_version=$(${CXX} -v 2>&1 | awk '/^.*clang version/ {print $3}')
			
			case ${clang_version} in
				4.*|5.*|6.*|7.*|8.*)
					AC_MSG_RESULT([yes])
					use_openmp="yes"
				;;
				*)
					AC_MSG_RESULT([no])
					AC_MSG_WARN([OpenMP unavailable])
					use_openmp="no"
				;;
			esac
		else
			AC_MSG_RESULT([no])
			AC_MSG_WARN([unsupported macOS build])
		fi
	fi
fi

openmp_flag=""
if test x"${use_openmp}" == x"yes"; then
	AC_MSG_CHECKING([for OpenMP])
	allldflags=$(${R_HOME}/bin/R CMD config --ldflags)
	hasOpenMP=$(echo ${allldflags} | grep -- -fopenmp)
	
	if test x"${hasOpenMP}" == x""; then
		AC_MSG_RESULT([no])
	else
		AC_MSG_RESULT([yes])
		openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
	fi
fi

AC_SUBST([OPENMP_CXXFLAGS], ["${openmp_flag}"])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
