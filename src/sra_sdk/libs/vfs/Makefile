# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================

default: all

TOP ?= $(shell ../../build/abspath.sh ../..)
MODULE = libs/vfs

INT_LIBS = \

EXT_LIBS = \
	libvfs

ALL_LIBS = \
	$(INT_LIBS) \
	$(EXT_LIBS)

include $(TOP)/build/Makefile.env

#-------------------------------------------------------------------------------
# outer targets
#
all std: makedirs
	@ $(MAKE_CMD) $(TARGDIR)/std

$(INT_LIBS): makedirs
	@ $(MAKE_CMD) $(ILIBDIR)/$@

$(EXT_LIBS): makedirs
	@ $(MAKE_CMD) $(LIBDIR)/$@

.PHONY: all std $(ALL_LIBS)

#-------------------------------------------------------------------------------
# std
#
$(TARGDIR)/std: \
	$(addprefix $(ILIBDIR)/,$(INT_LIBS)) \
	$(addprefix $(LIBDIR)/,$(EXT_LIBS))

.PHONY: $(TARGDIR)/std

#-------------------------------------------------------------------------------
# clean
#
clean: stdclean

.PHONY: clean

#-------------------------------------------------------------------------------
# tag
#
tag: \
	$(addsuffix _tag,$(ALL_LIBS))

.PHONY: tag $(addsuffix _tag,$(ALL_LIBS))

#-------------------------------------------------------------------------------
# vfs
#
$(LIBDIR)/libvfs: $(addprefix $(LIBDIR)/libvfs.,$(LIBEXT))

VFS_SRC = \
	path \
	syspath \
	manager 

VFS_OBJ = \
	$(addsuffix .$(LOBX),$(VFS_SRC))

VFS_LIB = \
	-dkrypto \
	-dkfg \
	-dkfs \
	-dklib

$(LIBDIR)/libvfs.$(SHLX): $(VFS_OBJ)
	$(LD) --dlib --vers $(SRCDIR) -o $@ $^ $(VFS_LIB)

$(LIBDIR)/libvfs.$(LIBX): $(VFS_OBJ)
	$(LD) --slib --vers $(SRCDIR) -o $@ $^ $(VFS_LIB)

libvfs_tag:
	@ $(TOP)/build/tag-module.sh $(MODULE) libvfs $(VFS_OBJ)



