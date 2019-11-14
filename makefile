Ejer.png : Ejer.dat Ejercicio28.py
	python Ejercicio28.py


Ejer.dat  : Ejer.x
	./Ejer.x 


Ejer.x : Ejercicio28.cpp
	c++ Ejercicio28.cpp -o Ejer.x
	

clean:
	rm Ejer.x Ejer.dat Ejer.png
