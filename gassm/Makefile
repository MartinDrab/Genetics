TARGET=gassm
BINDIR=../bin
CC=gcc
CFLAGS= -std=c99 -Wall
SHAREDDIR=../shared
TMPDIR=../tmp
OBJDIR=$(TMPDIR)/$(TARGET)
SHAREDOBJDIR=$(TMPDIR)/shared
LIB= -lm

OBJ=\
	$(OBJDIR)/gassm.o \
	$(OBJDIR)/test-file.o \
	$(OBJDIR)/assembly.o \
	$(OBJDIR)/kmer.o \
	$(OBJDIR)/kmer-edge.o \
	$(OBJDIR)/kmer-graph.o \
	$(OBJDIR)/kmer-table.o \
	$(SHAREDOBJDIR)/utils.o \
	$(SHAREDOBJDIR)/options.o

INCLUDE= -I./ -I./../include -I./$(SHAREDDIR)
SHAREDINCLUDE= -I./$(SHAREDDIR)

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

$(TARGET): $(OBJ)
	@mkdir -p $(BINDIR)
	@echo "Linking $@"
	@$(CC) $(OBJ) -o $(BINDIR)/$@ $(LIB)

all: $(TARGET)

clean:
	@rm -f $(OBJ)
	@rm -f $(SHAREDOBJ)