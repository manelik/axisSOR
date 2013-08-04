

vpath %.o objs
vpath %.f90 src

vpath sorxample exe

FF= gfortran

FILES := $(patsubst %.f90,%.o,$(wildcard src/*.f90))
FILES += $(patsubst %.f90,%.o,$(wildcard src/*/*.f90))

OBJS := $(notdir $(FILES))

%.o : %.f90
	@ echo "COMPILING FILE: $(notdir $<)"
	@ $(FF) -I objs -c $< -o objs/$@
	@ echo

.DEFAULT : ; @

start :  hello dir compile link 


hello :
	@ touch .timestart
	@ /bin/rm -f .config; echo $(ARCH) > .config
	@ echo
	@ echo "Compiling"
	@ echo
	@ echo

dir :	
	@ mkdir -p exe; mkdir -p objs

compile : $(OBJS)

link : sorxample

sorxample : $(OBJS)
	@ echo
	@ echo "LINKING ..."
	@ echo
	@ echo
	cd objs; $(FF) $(OBJS) -o ../exe/sorxample
	@ echo
	@ echo
	@ echo  "COMPILATION DONE!"
	@ echo
	@ echo

