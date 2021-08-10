# Introducción a la investigación Teórica 2021-1

Este repositorio contiene el codigo necesario para la investigación realizada
para la materia de introducción a la investigación teórica por:

-   Federico García
-   Ronald Cortés
-   Juan Bernardo Benavides

Estudiantes de física de la Universidad Nacional de Colombia sede Bogotá.

La simulación principal que se uso para la obtención de los resultados reportados
en el informe se encuentran en la carpeta `Codigo`. Estos se corrieron usando
computadores con sistema operativo Ubuntu, Mac con g++ y Python 3 instalados.

Para reproducir los resultados es necesario entrar a esta carpeta ` $ cd Codigo`
y correr el comando `make`. Este compilará la simulación en C++.

Este compilado funciona para sacar todos los resultados que se desee segun los
parametros adicionales que se escriban al correrlo. Para esto, en esa misma carpeta
hay multiples bash scripts que corren las simulaciones que se desean y organizan
los datos resultantes en una carpeta data.
El corazón del proyecto son los archivos funciones.cpp y main.cpp, el primero contiene
la información acerca de todas las funciones implementadas en la simulación, el segundo 
se compone de la función main con la implementación de dichas funciones, el ejecutable principal
es generado utilizando el comando `make main.x`. `main.x` recibe 6 parámetros cuyo orden son:
- 1000 kT
- 10 V
- 1000 R 
- $$-2/(\sigma-1)$$
- S
- Frontera
Se pueden correr simulaciones variando I, T y V (cada una es un parametro en
la simulación) y para cada uno existen diferentes variedades que se pueden correr.
Por ejemplo, existen 4 simulaciones diferentes (con diferentes parametros fijos)
en las que se obtiene la dimensión fractal en función de I. Para correr cada una
de estas se corre el make especificando la variable correspondiente al parametro.
Para I esto sería:

`$ make i=1 I1`

o

`$ make i=2 I2``

Es importante especificar la variable si es diferente de 1, ya que si no se especifica,
make dirá que el nombre no existe.

Después de esto, los datos quedan organizados en la carpeta data y en una carpeta
con el nombre correspondiente a la información que se corrió, i.e `Codigo/data/I1`.

Para graficar estos datos, se uso Python. Es necesario tener las siguientes librerías instaladas:

-   Numpy
-   Matplotlib
-   statsmodels

Se escoge cuales de las graficas se quiere sacar y esos datos se deben sacar de sus
carpetas individuales y pasarse directamente a la carpeta `data` que es donde Python
espera encontrarlos. Posteriormente en el archivo `Calcular.py` se descomentan las
lineas correspondientes a las graficas que se desea sacar y se corre `python Calcular.py`
o, si es necesario en su equipo `python3 Calcular.py`

Las graficas correspondientes estaran ahora en la carpeta `Codigo` en formato png.
Por ejemplo, si se desea realizar la gráfica de I2 basta con descomentar las siguientes lineas:
#generate_txt("imbalance",100, 1, 1)
#plot("Results_imbalance100T1V1R.txt","imbalance", "I2", "kT="+str(0.01)+", R/L=0.01, V=0.1")
Notemos que los argumentos de las funciones estan asociados con los parámetros recibidos por el main.x

