TARGET=gassm2
BINDIR=../bin
CC=gcc
CFLAGS= -Wall -std=c99
SHAREDDIR=../shared
LIBKMERDIR=../libkmer
TMPDIR=../tmp
OBJDIR=$(TMPDIR)/$(TARGET)
SHAREDOBJDIR=$(TMPDIR)/shared
LIBKMEROBJDIR=$(TMPDIR)/libkmer
LIBS= -lm

OBJ=\
	$(SHAREDOBJDIR)/utils.o \
	$(SHAREDOBJDIR)/options.o \
	$(SHAREDOBJDIR)/dym-array.o \
	$(SHAREDOBJDIR)/general-queue.o \
	$(SHAREDOBJDIR)/input-file.o \
	$(SHAREDOBJDIR)/paths.o \
	$(SHAREDOBJDIR)/reads.o \
	$(SHAREDOBJDIR)/ssw.o \
	$(OBJDIR)/gassm2.o \
	$(LIBKMEROBJDIR)/assembly.o \
	$(LIBKMEROBJDIR)/kmer-edge.o \
	$(LIBKMEROBJDIR)/kmer-graph.o \
	$(LIBKMEROBJDIR)/kmer-table.o \
	$(LIBKMEROBJDIR)/kmer.o \


INCLUDE= -I./ -I./../include -I./$(SHAREDDIR)

SHAREDINCLUDE= -I./../icnlude -I./$(SHAREDDIR)

$(OBJDIR)/%.o : %.c
	@mkdir -p $(TMPDIR)
	@mkdir -p $(OBJDIR)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(SHAREDOBJDIR)/%.o : $(SHAREDDIR)/%.c
	@mkdir -p $(TMPDIR)
	@mkdir -p $(SHAREDOBJDIR)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) $(SHAREDINCLUDE) -c $< -o $@

$(LIBKMEROBJDIR)/%.o : $(LIBKMERDIR)/%.c
	@mkdir -p $(TMPDIR)
	@mkdir -p $(SHAREDOBJDIR)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) $(SHAREDINCLUDE) -c $< -o $@

$(TARGET): $(OBJ)
	@mkdir -p $(BINDIR)
	@echo "Linking $@"
	@$(CC) $(OBJ) $(LIBS) -o $(BINDIR)/$@

all: $(TARGET)

clean:
	@rm -f $(OBJ)
	@rm -f $(SHAREDOBJ)