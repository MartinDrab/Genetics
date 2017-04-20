TARGET=vcffix
BINDIR=../bin
CC=g++
CFLAGS= -O3 -Wall --std=c++11 -D_FILE_OFFSET_BITS=64 -m64 -fopenmp -DNDEBUG -Wno-unused-function -g
TMPDIR=../tmp
OBJDIR=$(TMPDIR)/$(TARGET)
LIBS= -lm -fopenmp

OBJ=\
	$(OBJDIR)/main.o \

INCLUDE= -I./ -I./../include -Iquote

$(OBJDIR)/%.o : %.cpp
	@mkdir -p $(TMPDIR)
	@mkdir -p $(OBJDIR)
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