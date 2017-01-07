TARGET=gassm2
BINDIR=../bin
CC=gcc
CFLAGS= -O3 -Wall --std=gnu99 -D_FILE_OFFSET_BITS=64 -m64 -fopenmp -DNDEBUG -Wno-unused-function -g
SHAREDDIR=../shared
LIBKMERDIR=../libkmer
TMPDIR=../tmp
OBJDIR=$(TMPDIR)/$(TARGET)
SHAREDOBJDIR=$(TMPDIR)/shared
LIBKMEROBJDIR=$(TMPDIR)/libkmer
LIBS= -lm -fopenmp

OBJ=\
	$(SHAREDOBJDIR)/utils.o \
	$(SHAREDOBJDIR)/utils-lookaside.o \
	$(SHAREDOBJDIR)/file-utils.o \
	$(SHAREDOBJDIR)/options.o \
	$(SHAREDOBJDIR)/dym-array.o \
	$(SHAREDOBJDIR)/general-queue.o \
	$(SHAREDOBJDIR)/input-file.o \
	$(SHAREDOBJDIR)/reads.o \
	$(OBJDIR)/gassm2.o \
	$(LIBKMEROBJDIR)/ssw.o \
	$(LIBKMEROBJDIR)/paired-reads.o \
	$(LIBKMEROBJDIR)/assembly.o \
	$(LIBKMEROBJDIR)/found-sequence.o \
	$(LIBKMEROBJDIR)/read-info.o \
	$(LIBKMEROBJDIR)/kmer-edge.o \
	$(LIBKMEROBJDIR)/kmer-graph.o \
	$(LIBKMEROBJDIR)/kmer-table.o \
	$(LIBKMEROBJDIR)/kmer.o \
	$(LIBKMEROBJDIR)/variant-graph.o \



INCLUDE= -I./ -I./../include -I./$(SHAREDDIR) -I../libkmer -I-

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

$(LIBKMEROBJDIR)/%.o : $(LIBKMERDIR)/%.c
	@mkdir -p $(TMPDIR)
	@mkdir -p $(LIBKMEROBJDIR)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(TARGET): $(OBJ)
	@mkdir -p $(BINDIR)
	@echo "Linking $@"
	@$(CC) $(OBJ) $(LIBS) -o $(BINDIR)/$@

all: $(TARGET)

clean:
	@rm -f $(OBJ)
	@rm -f $(SHAREDOBJ)