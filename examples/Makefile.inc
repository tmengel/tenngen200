# TennGen configuration file.
# Generated on Wed 13 Apr 2022 11:51:35 PM EDT with the user supplied options:
# --cxx=g++
# --with-root-lib=/home/tmengel/root/lib
# --with-root-include=/home/tmengel/root/include

# Install directory prefixes.
PREFIX_BIN=/home/tmengel/tenngen200/bin
PREFIX_INCLUDE=/home/tmengel/tenngen200/include
PREFIX_LIB=/home/tmengel/tenngen200/lib
PREFIX_SHARE=/home/tmengel/tenngen200/share/TennGen200

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

FASTJET3_USE=false
FASTJET3_CONFIG=
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

ROOT_USE=true
ROOT_CONFIG=//home/tmengel/root/bin/root-config
ROOT_BIN=//home/tmengel/root/bin/
ROOT_INCLUDE=-I/home/tmengel/root/include
ROOT_LIB=-L/home/tmengel/root/lib -Wl,-rpath,/home/tmengel/root/lib -lCore
