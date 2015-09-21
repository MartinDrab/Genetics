TARGET=gassm
BINDIR=../bin
CC=gcc
CFLAGS= -Wall
SHAREDDIR=../shared
OBJDIR=../tmp/$(TARGET)
SHAREDOBJDIR=../tmp/shared

SOURCE=\
	gassm.c \
	kmer.c \
	kmer-edge.c \
	kmer-graph.c \
	kmer-table.c \

SHAREDSOURCE=\
	$(SHAREDDIR)\utils.c \
	$(SHAREDDIR)\options.c

OBJ=\
	$(OBJDIR)/gassm.o \
	$(OBJDIR)/kmer.o \
	$(OBJDIR)/kmer-edge.o \
	$(OBJDIR)/kmer-graph.o \
	$(OBJDIR)/kmer-table.o

SHAREDOBJ=\
	$(SHAREDOBJDIR)/utils.o \
	$(SHAREDOBJDIR)/options.o

INCLUDE= -I./ -I./../include -I./../shared

$(OBJDIR)/%.o : %.c
	@mkdir -p $(OBJDIR)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $(OBJDIR)/$@

$(SHAREDOBJDIR)/%.o : %.c
	@mkdir -p $(SHAREDOBJDIR)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) $(INCLUDE) -c $(SHAREDSOURCE)/$< -o $(SHAREDOBJDIR)/$@

$(TARGET): $(OBJ) $(SHAREDOBJ)
	@mkdir -p $(BINDIR)
	@echo "Linking $@"
	@$(CC) $(OBJ) $(SHAREDOBJ) -o $(BINDIR)/$@

all: $(TARGET)

clean:
	@rm -f $(OBJ)
	@rm -f $(SHAREDOBJ)