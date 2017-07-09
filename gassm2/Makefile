TARGET=gassm2
BINDIR=../bin
CC=gcc
CFLAGS= -O3 -Wall --std=gnu99 -D_FILE_OFFSET_BITS=64 -m64 -fopenmp -DNDEBUG -Wno-unused-function -g
SHAREDDIR=../shared
TMPDIR=../tmp
OBJDIR=$(TMPDIR)/$(TARGET)
SHAREDOBJDIR=$(TMPDIR)/shared
LIBKMEROBJDIR=$(TMPDIR)/libkmer
LIBS= -lm -fopenmp -L./../lib -l:librcorrect.a -l:libkmer.a

OBJ=\
	$(SHAREDOBJDIR)/utils.o \
	$(SHAREDOBJDIR)/utils-lookaside.o \
	$(SHAREDOBJDIR)/file-utils.o \
	$(SHAREDOBJDIR)/options.o \
	$(SHAREDOBJDIR)/input-file.o \
	$(SHAREDOBJDIR)/reads.o \
	$(SHAREDOBJDIR)/kthread.o \
	$(SHAREDOBJDIR)/ksw.o \
	$(OBJDIR)/gassm2.o \


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
	@$(CC) $^ $(LIBS) -o $(BINDIR)/$@

all: $(TARGET)

clean:
	@rm -f $(OBJ)
	@rm -f $(SHAREDOBJ)