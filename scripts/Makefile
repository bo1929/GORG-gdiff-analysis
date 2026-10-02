# Build helpers for the standalone C tools in this directory.
#
#   make -C scripts                   # build everything (OpenMP if available)
#   make -C scripts dip-pairs         # just the dip statistic tool
#   make -C scripts clean
#
# OpenMP is detected automatically; on macOS with Homebrew libomp the flags
# below are used, otherwise the tools build single-threaded (still correct,
# just serial).  Set OPENMP=0 to force a serial build.
#
# CFLAGS uses gnu99 rather than strict c99: the tools need POSIX getline(),
# ssize_t and off_t, which glibc hides under -std=c99.  dip_pairs.c also
# defines _POSIX_C_SOURCE/_DEFAULT_SOURCE itself, so -std=c99 works too, but
# gnu99 is the safer default on Linux.  Override CFLAGS to taste.

CC      ?= cc
CFLAGS  ?= -O3 -std=gnu99 -Wall -Wextra
LDLIBS  ?= -lm

OPENMP ?= auto
ifeq ($(OPENMP),auto)
  # probe: does the compiler accept -fopenmp as-is?
  OMP_OK := $(shell printf 'int main(void){return 0;}' > /tmp/.omp_probe_$$$$.c 2>/dev/null; \
                   $(CC) -fopenmp -o /tmp/.omp_probe_$$$$ /tmp/.omp_probe_$$$$.c >/dev/null 2>&1 \
                   && echo yes; rm -f /tmp/.omp_probe_$$$$*)
  ifeq ($(OMP_OK),yes)
    OMPFLAGS = -fopenmp
  else ifneq ($(wildcard /opt/homebrew/opt/libomp/include),)
    OMPFLAGS = -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
               -L/opt/homebrew/opt/libomp/lib -lomp
  else ifneq ($(wildcard /usr/local/opt/libomp/include),)
    OMPFLAGS = -Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include \
               -L/usr/local/opt/libomp/lib -lomp
  endif
else ifeq ($(OPENMP),0)
  OMPFLAGS =
else
  OMPFLAGS = -fopenmp
endif

all: dip-pairs dip-only estimate-sample-ani reconcile-gdiff-sample summarize-gdiff-sample

dip-pairs: dip_pairs.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -o $@ $< $(LDLIBS)

dip-only: dip_only.c dip_pairs.c
	$(CC) $(CFLAGS) -o $@ dip_only.c $(LDLIBS)

estimate-sample-ani: estimate_sample_ani.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

reconcile-gdiff-sample: reconcile-gdiff-sample.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

summarize-gdiff-sample: summarize-gdiff-sample.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f dip-pairs dip-only estimate-sample-ani reconcile-gdiff-sample summarize-gdiff-sample

.PHONY: all clean
