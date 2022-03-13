CC=g++
CFLAGS=-std=c++11 -O3 -fopenmp -g -mlzcnt
DLLFLAGS=-shared
LIBS=-llibpng16 -lgmp -lmpir -lfreeglut -lglew32 -lopengl32 -L.\libs\mpir -L.\libs\libpng
INC=-I.\png++ -I.\libs\mpir -I.\libs\libpng

all: core.dll mandel_gui mandel


core.dll: dll_interface.o dll_entry.o mandel_pertubation.o float_math.o
	$(CC) -o $@ $^ $(CFLAGS) $(DLLFLAGS) $(LIBS)

mandel_gui: gui_main.o mandel_navigator.o navigator.o glut_canvas.o mandel_pertubation.o image.o float_math.o gradient_window.o widget.o interpolate.o color.o observer.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

mandel: mandel.o mandel_pertubation.o image.o float_math.o renderer.o dll_interface.o interpolate.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)


mandel_pertubation.o: mandel_perturbation.cpp mandel_perturbation.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

mandel.o: mandel.cpp image.h mandel_perturbation.h test_cases.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

dll_interface.o: dll_interface.cpp dll_interface.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

renderer.o: renderer.cpp renderer.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

float_math.o: float_math.cpp float_math.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

image.o: image.cpp image.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

dll_entry.o: dll_entry.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

gui_main.o: gui_main.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

mandel_navigator.o: mandel_navigator.cpp mandel_navigator.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

navigator.o: navigator.cpp navigator.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

gradient_window.o: gradient_window.cpp gradient_window.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

interpolate.o: interpolate.cpp interpolate.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

widget.o: widget.cpp widget.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

color.o: color.cpp color.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

observer.o: observer.cpp observer.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

utils.o: utils.cpp utils.h
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)


.PHONY: clean

clean:
	rm -f *.o
