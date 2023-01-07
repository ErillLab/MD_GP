OBJS :=  ../src/_aux.o ../src/_multiplacement.o _organism.o _test.o

EXECUTABLE := test

CC := gcc

CFLAGS := -Wall -g

DEFS := 
INCLUDES :=
LDFLAGS := -lm

%.o: %.c %.h
	$(CC) $(CFLAGS) $(DEFS) $(INCLUDES) -c $< -o $@

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(DEFS) $(INCLUDES) $(OBJS) -o $@ $(LDFLAGS)
