TARGET=libkmer.a
BINDIR=../lib
CC=gcc
CFLAGS= -fpic -O3 -Wall --std=gnu99 -D_FILE_OFFSET_BITS=64 -m64 -fopenmp -DNDEBUG -Wno-unused-function -g
SHAREDDIR=../shared
TMPDIR=../tmp
OBJDIR=$(TMPDIR)/$(TARGET)
SHAREDOBJDIR=$(TMPDIR)/shared
LIBS= -lm -fopenmp

OBJ=\
	$(SHAREDOBJDIR)/utils.o \
	$(SHAREDOBJDIR)/utils-lookaside.o \
	$(SHAREDOBJDIR)/reads.o \
	$(SHAREDOBJDIR)/drand48.o \
	$(SHAREDOBJDIR)/kthread.o \
	$(OBJDIR)/assembly.o \
	$(OBJDIR)/variant.o \
	$(OBJDIR)/kmer-short.o \
	$(OBJDIR)/kmer-debug.o \
	$(OBJDIR)/kmer-edge.o \
	$(OBJDIR)/kmer-table.o \
	$(OBJDIR)/kmer-graph.o \
	$(OBJDIR)/read-info.o \
	$(OBJDIR)/paired-reads.o \
	$(OBJDIR)/ssw.o \
	$(OBJDIR)/variant-graph.o \

INCLUDE= -I./ -I./../include -I./$(SHAREDDIR) -I-

SHAREDINCLUDE= -I./../icnlude -I./$(SHAREDDIR) -I-

$(OBJDIR)/%.o : %.c
	@mkdir -p $(TMPDIR)
	@mkdir -p $(OBJDIR)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(SHAREDOBJDIR)/%.o : $(SHAREDDIR)/%.c
	@mkdir -p $(TMPDIR)
	@mkdir -p $(SHAREDOBJDIR)
	@echo "Compiling $<"
	@$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJ)
	@mkdir -p $(BINDIR)
	@echo "Linking $@"
	@ar rcs $(BINDIR)/$@ $^

all: $(TARGET)

clean:
	@rm -f $(OBJ)
	@rm -f $(SHAREDOBJ)