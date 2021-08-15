#########################################################################
#
# Makefile for cuda_test
#
#########################################################################

objects = \
adi-date-time.o \
adi-time-span.o \
adi-interval.o \
Util.o \
adi-type-define.o \
adi-util.o \
adi-constant.o \
adi-object.o \
adi-station.o \
adi-satellite.o \
adi-turntable.o \
adi-sun.o \
adi-link.o \
adi-station-kernel.o \
adi-satellite-kernel.o \
adi-sun-kernel.o \
adi-link-kernel.o \
adi-station-container.o \
adi-satellite-container.o \
adi-satellite-list.o \
adi-station-list.o \
adi-turntable-list.o \
adi-link-helper.o \
adi.o \
libadi.a \
adi-main.o
GCC_FLAG = -c -fPIC
NVCC_FLAG = --compiler-options "-Wall -Wfatal-errors -fPIC" -dc

edit:$(objects)

adi-date-time.o:adi-date-time.h adi-date-time.cc
	gcc $(GCC_FLAG) adi-date-time.cc

adi-time-span.o:adi-time-span.h adi-time-span.cc
	gcc $(GCC_FLAG) adi-time-span.cc

adi-interval.o:adi-interval.h adi-interval.cc
	gcc $(GCC_FLAG) adi-interval.cc

Util.o:Util.h Util.cc
	gcc $(GCC_FLAG) Util.cc

adi-util.o:adi-util.cu adi-util.h
	nvcc $(NVCC_FLAG) adi-util.cu

adi-constant.o:adi-constant.cc adi-constant.h
	nvcc $(NVCC_FLAG) adi-constant.cc

adi-type-define.o:adi-type-define.cc adi-type-define.h
	nvcc $(NVCC_FLAG) adi-type-define.cc
	
adi-object.o:adi-object.h adi-object.cu
	nvcc $(NVCC_FLAG) adi-object.cu

adi-station.o:adi-station.cu adi-station.h
	nvcc $(NVCC_FLAG) adi-station.cu

adi-satellite.o:adi-satellite.cu adi-satellite.h
	nvcc $(NVCC_FLAG) adi-satellite.cu

adi-turntable.o:adi-turntable.cc adi-turntable.h
	nvcc $(NVCC_FLAG) adi-turntable.cc

adi-sun.o:adi-sun.h adi-sun.cu
	nvcc $(NVCC_FLAG) adi-sun.cu

adi-link.o:adi-link.h adi-link.cu
	nvcc $(NVCC_FLAG) adi-link.cu

adi-station-kernel.o:adi-station-kernel.cu adi-station-kernel.h
	nvcc $(NVCC_FLAG) adi-station-kernel.cu

adi-satellite-kernel.o:adi-satellite-kernel.cu adi-satellite-kernel.h
	nvcc $(NVCC_FLAG) adi-satellite-kernel.cu
	
adi-sun-kernel.o:adi-sun-kernel.h adi-sun-kernel.cu
	nvcc $(NVCC_FLAG) adi-sun-kernel.cu

adi-link-kernel.o:adi-link-kernel.h adi-link-kernel.cu
	nvcc $(NVCC_FLAG) adi-link-kernel.cu

adi-satellite-list.o:adi-satellite-list.h adi-satellite-list.cc
	nvcc $(NVCC_FLAG) -c adi-satellite-list.cc

adi-station-list.o:adi-station-list.h adi-station-list.cc
	nvcc $(NVCC_FLAG) -c adi-station-list.cc

adi-turntable-list.o:adi-turntable-list.h adi-turntable-list.cc
	nvcc $(NVCC_FLAG)	-c adi-turntable-list.cc

adi-link-helper.o:adi-link-helper.h adi-link-helper.cu
	nvcc $(NVCC_FLAG) adi-link-helper.cu
	
adi-station-container.o:adi-station-container.h adi-station-container.cu
	nvcc $(NVCC_FLAG) adi-station-container.cu
	
adi-satellite-container.o:adi-satellite-container.h adi-satellite-container.cu
	nvcc $(NVCC_FLAG) adi-satellite-container.cu

adi.o: \
adi-date-time.o \
adi-time-span.o \
adi-interval.o \
Util.o \
adi-type-define.o \
adi-util.o \
adi-constant.o \
adi-object.o \
adi-station.o \
adi-satellite.o \
adi-turntable.o \
adi-sun.o \
adi-link.o \
adi-station-kernel.o \
adi-satellite-kernel.o \
adi-sun-kernel.o \
adi-link-kernel.o \
adi-station-container.o \
adi-satellite-container.o \
adi-satellite-list.o \
adi-station-list.o \
adi-turntable-list.o \
adi-link-helper.o
	nvcc --compiler-options "-Wall -Wfatal-errors -fPIC" -dlink -o $@ $^

libadi.a:adi.o
	ar -crv $@ *.o

adi-main.o:adi-main.cc libadi.a
	gcc -I/usr/local/cuda/include -o adi-main adi-main.cc -L. -L/usr/local/cuda/lib64 -ladi -lstdc++ -lm -lcudart -lcublas -lcurand

.PHONY:clean
clean:
	-rm -rf $(objects)