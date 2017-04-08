TARGET=rcorrect
BINDIR=../bin
CC=gcc
CFLAGS= -O3 -Wall --std=gnu99 -D_FILE_OFFSET_BITS=64 -m64 -fopenmp -DNDEBUG -Wno-unused-function -g
SHAREDDIR=../shared
TMPDIR=../tmp
OBJDIR=$(TMPDIR)/$(TARGET)
SHAREDOBJDIR=$(TMPDIR)/shared
LIBS= -lm -fopenmp

OBJ=\
	$(SHAREDOBJDIR)/utils.o \
	$(SHAREDOBJDIR)/file-utils.o \
	$(SHAREDOBJDIR)/input-file.o \
	$(SHAREDOBJDIR)/reads.o \
	$(OBJDIR)/bfc.o \
	$(OBJDIR)/bseq.o \
	$(OBJDIR)/htab.o \
	$(OBJDIR)/internal.o \
	$(OBJDIR)/main.o \


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