# Standalone build for bbe32ep_library (bypasses Xplorer).
#   make -f build.mk            build libnaturedsp_bbe32ep.a
#   make -f build.mk -j8        parallel
#   make -f build.mk clean

XT_TOOLS  := C:/usr/xtensa/XtDevTools/install/tools/RI-2020.5-win32/XtensaTools
XT_CORE   := bbe32_saturn
XT_SYSTEM := C:/usr/xtensa/XtDevTools/install/builds/RI-2020.5-win32/$(XT_CORE)/config

CC := $(XT_TOOLS)/bin/xt-clang
AR := $(XT_TOOLS)/bin/xt-ar

XTFLAGS := --xtensa-core=$(XT_CORE) --xtensa-system=$(XT_SYSTEM)
CFLAGS  := $(XTFLAGS) -std=c99 -O2 -Iinclude -Iinclude_private

OBJDIR := build/$(XT_CORE)
TARGET := $(OBJDIR)/libnaturedsp_bbe32ep.a

# Directories not built for bbe32_ep.
# fft/bcnfft is a hard exclude: 4 files reference BBE_SHFLI_MMC4X4X4X4_M2_STEP_2,
# a shuffle constant absent from this core's TIE headers, so they fail to compile.
LIB_EXCLUDE := fft/bcnfft matop/cmatvmul matop/cmatvmulf matop/cmatherm matop/cmathermf \
               matop/rcmatmul matop/rcmatmulf matop/rcmatvmul matop/rcmatvmulf \
               matop/matmul matop/matmulf matop/cmatmul matop/cmatmulf \
               matinv/bcholn matinv/bcholnf

SRC_DIRS := comm complex complexf fft fir fit id iir iirf math mathf \
            matinv matop tables vector vectorf

# Use make's own wildcard rather than shell find: on Windows the shell make
# spawns strips the quotes around '*.c', so it globs in the CWD before find
# runs. Sources sit at most 2 levels below a top dir (e.g. fft/bcfft/x.c).
ALL_SRCS := $(wildcard feature.c version.c) \
            $(foreach d,$(SRC_DIRS),$(wildcard $(d)/*.c $(d)/*/*.c $(d)/*/*/*.c))
# Anchor each pattern at a path separator so matop/matmul does not also
# swallow matop/matmulf; the trailing / makes the prefix unambiguous.
SRCS := $(filter-out $(foreach d,$(LIB_EXCLUDE),$(d)/%),$(ALL_SRCS))
OBJS := $(addprefix $(OBJDIR)/,$(SRCS:.c=.o))

.PHONY: all clean list
all: $(TARGET)

$(TARGET): $(OBJS)
	@echo "AR  $@ ($(words $(OBJS)) objects)"
	@rm -f $@
	@find $(OBJDIR) -name '*.o' | xargs -n 200 $(AR) rcs $@

$(OBJDIR)/%.o: %.c
	@mkdir -p $(dir $@)
	@$(CC) $(CFLAGS) -c $< -o $@

list:
	@echo "excluded dirs : $(words $(LIB_EXCLUDE))"
	@echo "sources total : $(words $(ALL_SRCS))"
	@echo "sources built : $(words $(SRCS))"
	@echo "skipped       : $(words $(filter-out $(SRCS),$(ALL_SRCS)))"

clean:
	rm -rf $(OBJDIR)

check:
	@echo "matmul kept   : $(words $(filter matop/matmul/%,$(SRCS)))  (want 0)"
	@echo "matmulf kept  : $(words $(filter matop/matmulf/%,$(SRCS)))  (want 0)"
	@echo "cmatvmul kept : $(words $(filter matop/cmatvmul/%,$(SRCS)))  (want 0)"
	@echo "bcfft kept    : $(words $(filter fft/bcfft/%,$(SRCS)))  (want 8, not excluded)"
	@echo "bcnfft kept   : $(words $(filter fft/bcnfft/%,$(SRCS)))  (want 0)"
	@echo "cholnf kept   : $(words $(filter matinv/cholnf/%,$(SRCS)))  (want >0, only bcholnf excluded)"
