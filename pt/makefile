PROJECT = pt

SRCS = src/main.c

DEBUG_FLAGS = -g \
              -fsanitize=address \
              -fsanitize=leak \
              -fsanitize=undefined \
              -fsanitize=bool \
              -fsanitize=enum \
              -fsanitize=float-cast-overflow \
              -fsanitize=signed-integer-overflow

RELEASE_FLAGS = -O3

$(PROJECT): $(SRCS)
	$(CC) -o $(PROJECT) $(SRCS) $(RELEASE_FLAGS) -L../ -I../ -l:sbmf.a -lm -lpthread -fopenmp -lgfortran
