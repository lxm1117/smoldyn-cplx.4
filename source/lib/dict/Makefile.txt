CC=gcc
CFLAGS=-g3 -Wall

all: test_dict

test: all
	./test_dict

test_dict: test_dict.o dict.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	$(RM) test_dict *.o
