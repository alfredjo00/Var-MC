CC = gcc

CFLAGS = \
	 -g \
	 -Iinclude \
	 -Wall \
	 -Werror \
	 -pedantic \
	 -fsanitize=address \
	 -fno-omit-frame-pointer

CFLAGS_OPT = \
	     -O2 \
	     -march=native

LIB = \
      -lm \
      -lgsl \
      -lgslcblas

OBJ = \
      obj/tools.o \
      obj/run.o  \
      obj/metropolis.o
	  
MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

H2: obj _H2

_H2: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o H2 $(LIB)

obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

obj:
	mkdir -p obj

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f program run-test
	rm -rf obj obj_test
	rm -f H2

.PHONY: clean
