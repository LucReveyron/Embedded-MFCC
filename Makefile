# Compilers
CXX = g++
CC = gcc

# Flags
CXXFLAGS = -Wall -g -std=c++17 -DKISSFFT_DATATYPE=int16_t -I./kissFFT
CFLAGS   = -Wall -g -DKISSFFT_DATATYPE=int16_t -I./kissFFT

# Target
TARGET = main

# Sources
CXX_SRCS = main.cpp
C_SRCS  = kissFFT/kiss_fft.c kissFFT/kiss_fftr.c

# Objects
CXX_OBJS = $(CXX_SRCS:.cpp=.o)
C_OBJS   = $(C_SRCS:.c=.o)
OBJS     = $(CXX_OBJS) $(C_OBJS)

# Default
all: $(TARGET) run

# Compile C++ files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile C files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Link
$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS)

# Run
run: $(TARGET)
	./$(TARGET)

# Clean
clean:
	rm -f $(TARGET) $(OBJS)
