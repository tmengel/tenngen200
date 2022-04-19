################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script
# and the distribution structure.
################################################################################

# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration and set the local directory structure.
ifeq (,$(findstring clean, $(MAKECMDGOALS)))
  ifneq (Makefile.inc,$(MAKECMDGOALS))
    LOCAL_CFG:=$(shell $(MAKE) Makefile.inc))
    include Makefile.inc
  endif
endif
LOCAL_BIN=bin/
LOCAL_DOCS=README.md ../../examples/Makefile.inc
LOCAL_EXAMPLE=examples
LOCAL_INCLUDE=include
LOCAL_LIB=lib
LOCAL_SHARE=share/TennGen200
LOCAL_SRC=src
LOCAL_TMP=tmp
LOCAL_MKDIRS:=$(shell mkdir -p $(LOCAL_TMP) $(LOCAL_LIB))
CXX_COMMON:=-I$(LOCAL_INCLUDE) $(CXX_COMMON)
OBJ_COMMON:=-MD $(CXX_COMMON) $(OBJ_COMMON) $(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`
LIB_COMMON=-Wl,-rpath,../lib:$(PREFIX_LIB) -ldl $(GZIP_LIB)

# TENNGEN.

OBJECTS=$(patsubst $(LOCAL_SRC)/%.cpp,$(LOCAL_TMP)/%.o,\
	$(sort $(wildcard $(LOCAL_SRC)/*.cpp)))
TARGETS= $(LOCAL_LIB)/libtenngen200.a $(LOCAL_LIB)/libtenngen200$(LIB_SUFFIX) #$(LOCAL_TMP)/TennGen200Data.o  

################################################################################
# RULES: Definition of the rules used to build TENNGEN.
################################################################################

# Rules without physical targets (secondary expansion for documentation).
.SECONDEXPANSION:
.PHONY: all install distclean

# All targets.
all: $(TARGETS) $(addprefix $(LOCAL_SHARE)/, $(LOCAL_DOCS)) $(LOCAL_SHARE)/AuAu200Data.root

# The documentation.
$(addprefix $(LOCAL_SHARE)/, $(LOCAL_DOCS)): $$(notdir $$@)
	cp $^ $@

# The Makefile configuration.
Makefile.inc:
	./configure


# TENNGEN.

$(LOCAL_TMP)/TennGen.o: $(LOCAL_SRC)/TennGen.cpp Makefile.inc
	$(CXX) $< -o $@ -c $(OBJ_COMMON) 
$(LOCAL_TMP)/%.o: $(LOCAL_SRC)/%.cpp
	$(CXX) $< -o $@ -c $(OBJ_COMMON)
$(LOCAL_LIB)/libtenngen200.a: $(OBJECTS)
	ar cr $@ $^
$(LOCAL_LIB)/libtenngen200$(LIB_SUFFIX): $(OBJECTS)
	$(CXX) $^ -o $@ $(CXX_COMMON) $(CXX_SHARED) $(CXX_SONAME)$(notdir $@)\
	  $(LIB_COMMON)

CXX_COMMON:=-I$(LOCAL_INCLUDE) $(CXX_COMMON) $(GZIP_LIB)
CXX_COMMON+= -L$(LOCAL_LIB) -Wl,-rpath,$(LOCAL_LIB) -ltenngen200 -ldl
TENNGEN=$(LOCAL_LIB)/libtenngen200$(LIB_SUFFIX)
$(LOCAL_SHARE)/AuAu200Data.root: $(TENNGEN) $(LOCAL_SRC)/TennGen200Data.cpp
	$(CXX) $< -o AuAuData $(LOCAL_SRC)/TennGen200Data.cpp -w $(CXX_COMMON) $(ROOT_LIB)\
	 `$(ROOT_CONFIG) --cflags --glibs`
	./AuAuData AuAu200Data.root 
	cp AuAu200Data.root $(LOCAL_SHARE)/AuAu200Data.root
	cp -rf distro-PNG $(LOCAL_SHARE)/
	cp AuAuData $(LOCAL_SHARE)/AuAuData
	rm -rf AuAu200Data.root AuAuData distro-PNG

# Install.
install: all
	mkdir -p $(PREFIX_BIN) $(PREFIX_INCLUDE) $(PREFIX_LIB) $(PREFIX_SHARE)
	rsync -a $(LOCAL_BIN)/* $(PREFIX_BIN)
	rsync -a $(LOCAL_INCLUDE)/* $(PREFIX_INCLUDE)
	rsync -a $(LOCAL_LIB)/* $(PREFIX_LIB)
	rsync -a $(LOCAL_SHARE)/* $(PREFIX_SHARE)
	rsync -a $(LOCAL_EXAMPLE) $(PREFIX_SHARE)

# Clean.
clean:
	rm -rf $(LOCAL_TMP) $(LOCAL_LIB)
	rm -f $(LOCAL_EXAMPLE)/*.root
	rm -f $(LOCAL_EXAMPLE)/tenngen-*
	rm -f $(LOCAL_EXAMPLE)/*.inc

# Clean all temporary and generated files.
distclean: clean
	cd examples && make clean
	find . -type f -name Makefile.inc -print0 | xargs -0 rm -f
	rm -rf $(LOCAL_BIN)
	rm -f $(LOCAL_SHARE)/README.md
	rm -f $(LOCAL_SHARE)/AuAuData
	rm -f $(LOCAL_SHARE)/AuAu200Data.root
