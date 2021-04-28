################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/match_gap_coro.cpp 

OBJS += \
./src/match_gap_coro.o 

CPP_DEPS += \
./src/match_gap_coro.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	clang++ -std=c++20 -I. -O3 -g3 -Wall -c -mllvm -inline-threshold=16384 -fmessage-length=0 -stdlib=libc++  -fcoroutines-ts -fpermissive -fno-exceptions -pthread -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


