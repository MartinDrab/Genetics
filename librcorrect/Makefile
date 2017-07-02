TARGET=librcorrect.a
BINDIR=../lib
CC=gcc
CFLAGS= -O3 -Wall --std=gnu99 -D_FILE_OFFSET_BITS=64 -m64 -fopenmp -DNDEBUG -Wno-unused-function -g
SHAREDDIR=../shared
TMPDIR=../tmp
OBJDIR=$(TMPDIR)/$(TARGET)
SHAREDOBJDIR=$(TMPDIR)/shared
LIBS= -lm -fopenmp

OBJ=\
	$(SHAREDOBJDIR)/utils.o \
	$(SHAREDOBJDIR)/reads.o \
	$(SHAREDOBJDIR)/drand48.o \
	$(SHAREDOBJDIR)/kthread.o \
	$(OBJDIR)/bfc.o \
	$(OBJDIR)/bseq.o \
	$(OBJDIR)/htab.o \
	$(OBJDIR)/internal.o \
	$(OBJDIR)/librcorrect.o \


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
	@ar rcs -o $(BINDIR)/$@ $<

all: $(TARGET)

clean:
	@rm -f $(OBJ)
	@rm -f $(SHAREDOBJ)