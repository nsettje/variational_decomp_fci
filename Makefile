#
# Plugin Makefile generated by Psi4.
#
# You shouldn't need to modify anything in this file.
#

# Location of your PSI4 source
top_srcdir = /home/settje/psi4
# Location of your PSI4 install, by default as listed
top_objdir = /home/settje/psi4/obj

# Start by figuring out whether we're on Linux or Mac (sorry, Mr. Gates)
UNAME := $(shell uname)
include $(top_objdir)/src/bin/MakeVars

# Reset these values, MakeVars changes them to valud only valid in Psi4's objdir
# Location of your PSI4 source
top_srcdir = /home/settje/psi4
# Location of your PSI4 install, by default as listed
top_objdir = /home/settje/psi4/obj

PSITARGET = $(shell basename `pwd`).so
PSILIBS = -L$(top_objdir)/lib -lPSI_plugin

CXXSRC = $(notdir $(wildcard *.cc))
#FXXSRC = $(notdir $(wildcard *.F))
DEPENDINCLUDE = $(notdir $(wildcard *.h))

BINOBJ = $(CXXSRC:%.cc=%.o) nosd_bfgs.o
default:: $(PSITARGET)

# Add the flags needed for shared library creation
ifeq ($(UNAME), Linux)
    LDFLAGS = -shared -fopenmp
endif
ifeq ($(UNAME), Darwin)
    LDFLAGS = -shared -undefined dynamic_lookup
    CXXOTH += -fno-common
endif
# The object files
%.o: %.cc
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -fPIC -c $< $(OUTPUT_OPTION)
%.o: %.F
	$(FC) -c -fPIC $< $(OUTPUT_OPTION)

$(PSITARGET): $(BINOBJ)
	$(CXX) $(LDFLAGS) -fPIC -o $@ $^ $(CXXDBG) $(PSILIBS)

	@echo "$(BINOBJ)"
# Erase all compiled intermediate files
clean:
	rm -f $(BINOBJ) $(PSITARGET) *.d *.pyc

# Dependency handling
%.d: %.cc
	$(CXXDEPEND) $(CXXDEPENDFLAGS) $(CXXFLAGS) $(CXXINCLUDE) $< | sed 's/^$*.o/$*.o $*.d/g' > $(@F)

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

