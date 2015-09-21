TARGET=gengen
BINDIR=../bin
CC=gcc
CFLAGS= -Wall -std=c99
SHAREDDIR=../shared
TMPDIR=../tmp
OBJDIR=..$(TMPDIR)/$(TARGET)
SHAREDOBJDIR=$(TMPDIR)/shared
LIBS= -lm

OBJ=\
	$(SHAREDOBJDIR)/utils.o \
	$(SHAREDOBJDIR)/options.o \
	$(OBJDIR)/generator.o \
	$(OBJDIR)/gengen.o \


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

$(TARGET): $(OBJ)
	@mkdir -p $(BINDIR)
	@echo "Linking $@"
	@$(CC) $(OBJ) $(LIBS) -o $(BINDIR)/$@

all: $(TARGET)

clean:
	@rm -f $(OBJ)
	@rm -f $(SHAREDOBJ)