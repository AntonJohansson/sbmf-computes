PROJECT = gp

SRCS = src/main.c

DEBUG_FLAGS = -g \
              -fsanitize=address \
              -fsanitize=leak \
              -fsanitize=undefined \
              -fsanitize=bool \
              -fsanitize=enum \
              -fsanitize=float-cast-overflow \
              -fsanitize=signed-integer-overflow

$(PROJECT): $(SRCS)
	$(CC) -o $(PROJECT) $(SRCS) -O3 -L../ -I../ -l:sbmf.a -lm -lpthread -fopenmp -lgfortran -L/home/aj/.local/lib/ -I/home/aj/.local/include
