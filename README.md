El archivo oz.py genera la solución de 0 a 20 partículas para todos los valores del parámetro s
Al utilizar range en los ciclos, es necesario poner el valor final al que queremos llegar +1 

Estas soluciones se cargan en ggrande2.py  para obtener las soluciones para todos los números de partículas para todos los valores del parámetro s
Los archivos se guardan en un archivo nombrado de la siguiente forma: 's*100Ngmcgb.txt' Por ejemplo para s=0.95 y 100 partículas, el archivo es 95100gmcgb.txt
La primera columa es el valor de gamma(theta), la segunda el de c(theta), lueego el de g(theta) y lluego el de b(theta). También se generan las soluciones para un punto 
anterior y un punto después las cuales ocuparemos para el calculo de la derivada numérica al calcular la consistencia termodinámica.

El archivo graficacgb.py nos hace las gráficas para gamma(theta), c(theta) y g(theta) para todos los valores de s y N y los gaurda en sus correspondientes carpetas

En el archivo consistencia_termo.py se cargan las soluciones para g(theta) obtenidas y se realizan los correspondientes cálculos de las ecuaciones de la consistencia termodnámica

El archivo grafica.py genera las gráficas de los puntos de cruce.
