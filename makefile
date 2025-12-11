obj = main.cpp lib.cpp lib.h args.h mytime.cpp mytime.h
exec = a.out ga.out gp.out gpO3.out tsan
cxx = -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
parallel = -pthread  $(cxx)
<<<<<<< HEAD
tsan = -fsanitize=thread $(parallel)
=======
tsan =  $(parallel) -fsanitize=thread
>>>>>>> 813702e ('Thu Dec 11 23:11:35 ')
gprO3 = -pg $(cxx)
gpr = -pg $(wo3)
wo3 = -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
testfold =matrixtests/
testexe = ./a.out



a.out: $(obj)

	g++ $(parallel) $(obj) -o a.out

tsan: $(obj)

	g++ $(tsan) $(obj) -o $@

mainzip:
	zip Frolov_PS.zip $(obj) makefile 

zip:
	zip pll.zip $(obj) makefile  


gp.out: $(obj)

	g++ $(gpr) $(obj) -o gp.out

gpO3.out: $(obj)

	g++ $(gprO3) $(obj) -o gpO3.out

ga.out: $(obj)

	g++ $(wo3) $(obj) -o ga.out

clean:
	rm $(exec)

git:
	../sequentialGauss/commit.sh
	
	

push:
	git push

test:
	./test.sh  $(testexe) $(testfold)

all: $(exec)
