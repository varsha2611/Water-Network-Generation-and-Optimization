UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	#echo $(UNAME)
	INTER = python   
	CMDL  = setup.py build_ext --inplace
else
	#XTRA_INCLUDES=C:\Py25\Lib\site-packages\numpy\core\include\numpy
	#export INCLUDE=$(INCLUDES);$(XTRA_INCLUDES)
	INTER = python   
	CMDL  = setup.py build_ext --inplace -c mingw32
endif

all: resMetricX

#resMetricX: resMetricX.pyx
resMetricX:
	#$(INTER) $(CMDL) $(OPTNS)
	$(INTER) $(CMDL)
clean:
	rm -f *.c
	rm -f *.pyd
	rm -f *.so
	rm -rf build

