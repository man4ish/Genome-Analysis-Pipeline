# $Header: /Synamatix/Projects/Repository/synabase/sxapp/src/Mkgcc.mak,v 1.45 2008/03/26 08:47:45 theanhoo Exp $

# OSTYPE variable check.
RESULTUNKNOWN=1
ifeq ($(OSTYPE),linux-gnu)
   RESULTUNKNOWN=0
endif
ifeq ($(OSTYPE),cygwin)
   RESULTUNKNOWN=0
endif
ifeq ($(OSTYPE),msys)
   RESULTUNKNOWN=0
endif
ifeq ($(RESULTUNKNOWN),1)
   ifeq ($(shell uname),Linux)
      OSTYPE=linux-gnu
      RESULTUNKNOWN=0
   endif
endif
ifeq ($(RESULTUNKNOWN),1)
   ifeq ($(TERM),cygwin)
      OSTYPE=cygwin
      RESULTUNKNOWN=0
   endif
endif
ifeq ($(RESULTUNKNOWN),1)
   $(error OSTYPE '$(OSTYPE)' not known)
endif

# Declare file suffix for gcc -E output files.
.SUFFIXES: .e

# Declare command aliases.
ECHO=echo
CP=cp

# Declare 3rd party API include & library paths.
REFPATH=/opt/synamatix/synabase/synabase
XERCESCLIB=$(REFPATH)/../xercesc/lib
SXAPIINCDIR=$(REFPATH)/include
SXAPILIBDIR=$(REFPATH)/lib

# Declare general & local compile flags.
CLFLAGS=-pedantic -Wall -Wno-deprecated -D_REENTRANT -D_LARGEFILE_SOURCE=1 -D_USE_FILE_OFFSET=64 -fPIC
ifeq ($(OSTYPE),linux-gnu)
   ifeq ($(HOSTTYPE),ia64)
      CLFLAGS=-Wno-deprecated -D_REENTRANT -D_LARGEFILE_SOURCE=1 -D_USE_FILE_OFFSET=64 -fPIC
   else
      CLFLAGS=-Wno-deprecated -D_REENTRANT -fPIC
   endif
else
   ifeq ($(OSTYPE),cygwin)
      CLFLAGS=-Wno-deprecated -D_REENTRANT
   endif
   ifeq ($(OSTYPE),msys)
      CLFLAGS=-Wno-deprecated -D_REENTRANT
   endif
endif

# Declare externally overridable include paths.
#COINC=
# Declare platform/compiler/machine specific include paths.
CSINC=
# Declare general & local include paths.
CLINC=-I. -I$(SXAPIINCDIR) 
# Declare general & local link flags.
BINRUNPATH=-Wl,-rpath,'$$LD_RUN_PATH',-rpath,'$$SYNABASEROOT/lib',-rpath,'$$ORIGIN/lib',-rpath,.,-rpath,/usr/local/lib
LDLFLAGS=$(BINRUNPATH)
# Declare externally overridable library paths.
#LDOINC=
# Declare platform/compiler/machine specific library paths.
LDSINC=
# Declare general & local library paths.
LDLINC=-L$(SXAPILIBDIR) -L$(XERCESCLIB)

# Declare executable file output path.
#XAPPPATH=/home/syhwah/bin
# Declare object file output path.
OBJOUTPATH=../obj
# Declare dependent libraries.
ifeq ($(findstring icc,$(CC)),icc)
LIBSTDCPP=-lcxa
LIBS=-i-static -static-libcxa -lsynabase -lneurabase -lxerces-c -lsxutilities $(LIBSTDCPP) -lm
else
LIBSTDCPP=-lstdc++
LIBS= -lsynabase -lneurabase -lxerces-c -lsxutilities $(LIBSTDCPP) -lm
endif

.c.e:
	-@$(ECHO) "$@: "
	$(CC) $(CFLAGS) $(COFLAGS) $(CLFLAGS) $(COINC) $(CSINC) $(CLINC) -dD -E -c $? > $@

.c.o:
	-@$(ECHO) "$@: "
	$(CC) $(CFLAGS) $(COFLAGS) $(CLFLAGS) $(COINC) $(CSINC) $(CLINC) -o $@ -c $?

$(OBJOUTPATH)/%.o: %.c
	-@$(ECHO) "$@: "
	$(CC) $(CFLAGS) $(COFLAGS) $(CLFLAGS) $(COINC) $(CSINC) $(CLINC) -o $@ -c $?

.cpp.e:
	-@$(ECHO) "$@: "
	$(CXX) $(CFLAGS) $(COFLAGS) $(CLFLAGS) $(COINC) $(CSINC) $(CLINC) -dD -E -c $? > $@ 

.cpp.o:
	-@$(ECHO) "$@: "
	$(CXX) $(CFLAGS) $(COFLAGS) $(CLFLAGS) $(COINC) $(CSINC) $(CLINC) -o $@ -c $?

$(OBJOUTPATH)/%.o: %.cpp
	-@$(ECHO) "$@: "
	$(CXX) $(CFLAGS) $(COFLAGS) $(CLFLAGS) $(COINC) $(CSINC) $(CLINC) -o $@ -c $?

default: all

all:
	-@$(ECHO) "$@: "
	$(CXX) $(CFLAGS) $(COFLAGS) $(CLFLAGS) $(COINC) $(CSINC) $(CLINC) $(LDFLAGS) $(LDOFLAGS) $(LDLFLAGS) $(LDOINC) $(LDSINC) $(LDLINC) -o SXSynoMap SXSynoMap.cpp $(LIBS)


clean:
	-@$(ECHO) "$@: "
	-@rm -vRf $(OBJOUTPATH)/*.o
	-@rm *.~*
	-@rm *.$*
	-@rm *.tmp
	-@rm *.bak
	-@rm *.tds
	-@rm *.e
	-@rm *.obj
	-@rm *.o
	-@rm *.lib
	-@rm *.a
	-@rm *.dll
	-@rm *.so*
	-@rm *.exe

veryclean: clean
	-@$(ECHO) "$@: "

# Generate release build (dynamically linked)
reld: COFLAGS=-O3 -DNDEBUG
reld: LDOFLAGS=-Wl,-S,-x
reld: all
	-@$(ECHO) "$@: "

# Generate release build (statically linked)
rels: COFLAGS=-O3 -DNDEBUG
rels: LDOFLAGS=-Wl,-S,-x
ifeq ($(OSTYPE),linux-gnu)
rels: LIBS=-lsxaddons -lsynabase -lxerces-c -lsxutilities -Wl,-Bstatic,$(LIBSTDCPP),-Bdynamic -lm
endif
rels: all
	-@$(ECHO) "$@: "

# Generate development/debug build (dynamically linked)
devd: COFLAGS=-g
devd: all
	-@$(ECHO) "$@: "

# Generate development/debug build (statically linked)
devs: COFLAGS=-g
ifeq ($(OSTYPE),linux-gnu)
devs: LIBS=-lsxaddons -lsynabase -lxerces-c -lsxutilities -Wl,-Bstatic,$(LIBSTDCPP),-Bdynamic -lm
endif
devs: all
	-@$(ECHO) "$@: "
