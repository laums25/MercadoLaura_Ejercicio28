Ejercicio28.png : Ejercicio28.dat Ejercicio28.py
	python Ejercicio28.py

Ejercicio28.dat  : Ejer.x
	./Ejer.x 

Ejer.x : Ejercicio28.cpp
	c++ Ejercicio28.cpp -o Ejer.x
	
clean:
	rm Ejer.x Ejercicio28.dat Ejercicio28.png
