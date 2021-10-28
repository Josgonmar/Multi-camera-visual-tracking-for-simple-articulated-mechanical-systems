## Multi-camera-visual-tracking
# Español:
Este proyecto tiene como objetivo asentar las bases para el desarrollo de un método sencillo y de bajo coste, capaz de localizar tanto en posición como en orientación un sistema mecánico articulado mediante la aplicación de técnicas de localización 3D basadas en un sistema multicámara calibrado, que consiste en un cubículo de 3x3x3 metros y cuenta con cuatro cámaras infrarrojas, en posiciones enfrentadas dos a dos.
El algoritmo, desarrollado en lenguaje C++ para Visual Studio, permite mediante la aplicación de técnicas de visión por computador como la Triangulación y el problema de la Perspectiva de N Puntos (PnP) localizar una serie de marcadores pasivos en el espacio tridimensional. Una distribución de dichos marcadores en puntos clave, junto con un modelo cinemático conocido de la estructura, permite estimar todos los grados de libertad de la misma aplicando mínimos cuadrados sobre una función que relaciona la localización tridimensional de los marcadores con el modelo cinemático. Se muestran resultados experimentales para un conjunto de estructuras de prueba con diferente número de grados de libertad, que permiten evaluar la eficacia del método propuesto y sus posibles mejoras para trabajos futuros.
Este problema tiene aplicación para el seguimiento del movimiendo de sistemas robóticos y vehículos articulados, así como para la captura de movimiento del cuerpo humano o partes de él.

FUNCIONAMIENTO:
La idea es capturar 3 imágenes desde puntos diferentes, de cierta articulación móvil cuya posición relativa a las cámaras es desconocida.
Este brazo articulado debe poseer marcadores (en este caso, de tamaños diferenciables) en cada una de las uniones entre eslabones, asi como en la que sería la base del mismo.
De esta forma, se detectan los centros de los marcadores en cada una de las tres imágenes aplicando una serie de operaciones morfológicas simples, y mediante triangulación se obtiene la posición global de cada uno de estos marcadores, respecto de una de las cámaras que se eligirá como referencia. Por supuesto, se entiende que dicho sistema de cámaras se encuentra completamente calibrado. Con las coordenadas trianguladas, se ha de resolver mediante mínimos cuadrados, una ecuación que relaciona dichas coordenadas con las coordenadas relativas de cada marcador, vistas desde la base del brazo. Esta ecuación, o sistema de ecuaciones, no es más que una transformación homogénea de coordenadas, por lo que obtenemos un total de 3 ecuaciones por cada marcador, y (3+X) incógnitas, donde X corresponde con el número de grados de libertad del brazo articulado, quedando de la sigueinte manera: (Xc,Yc,Zc)=F(R,P,Y,Qn), teniendo en cuenta que también nos interesa conocer los angulos de giro a aplicar para obtener tanto la posición como la orientación relativa respecto de la cámara central.
La imágen resultante muestra el valor de los ángulos de alabeo, cabeceo y guiñada, así como de los grados de libertad, ademas de una representación gráfica de como quedarían dispuestos los ejes de referencia del brazo articulado.

COSAS A TENER EN CUENTA:
- Por motivos ajenos, no fue posible trabajar con las 4 cámaras a la vez, por lo que está pensado para funcionar con solo 3 de ellas. Más camaras incurre en, obviamente, una mejora en la presición de las estimaciones.
- La función que opera la triangulación contenida en la librería OpenCV solo es capaz de trabajar con 2 cámaras a la vez, por lo que habrá que promediar los resultados entre todas.
- Los progrmas están pensados para funcionar con una configuración muy concreta de los sistemas articulados, de 2, 3 y 4 grados de libertad, con movimientos en 3 dimensiones. Se supone que tanto las dimensiones como las ecuaciones cinemáticas son conocidas. No obstante, el sistema esta pensado para fucionar con cualquier configuración, siempre y cuando se modifiquen las partes correspondientes.
- Si bien se puede restringir matematicamente a que la resolución siempre encuentre uno de los ejes de referencia como solidario a uno de los eslabones, no es posible hacerlo con los restantes, por lo que no existe una solución única al problema. Pero, normalmente el error de reproyección no suele ser superior a 3 píxeles.

REQUISITOS:

- El código está escrito íntegramente en lenguaje C/C++.
- OpenCV 4.5.2
- Eigen 3.3.9
- Ceres 2.0

# English:
This project aims to settle the foundations for the development of a simple and low cost method, capable of locating both, position and orientation, of an articulated mechanical system through the application of 3D localization techniques based on a calibrated multi camera system, consisting of a cubicle of 3x3x3 meters and four infrared cameras, placed in opposing positions two by two.
The algorithm, developed in C++ language for Visual Studio, allows by applying some computer vision techniques such as Triangulation and the Perspective N Point (PnP) problem to locate a series of passive markers in a three dimensional space. A distribution of these markers at key points, together with a known kinematic model of the structure, enables to estimate all the degrees of freedom by applying least squares on a function that relates the three dimensional location of the markers with the kinematic model. Experimental results are shown for a set of test structures with different number of degrees of freedom, which allow to assess the effectiveness of the proposed method and its possible improvements for future work.
This problem has application for the tracking of the movement of robotic systems and artificial vehicles, as well as for the capture of movement of the human body or parts of it.

HOW IT WORKS:
The idea its to capture three images of a certain articulated mechanical system, whose position and orientation from the cameras is unknown.
This articulated arm is supposed to have several markers (and to be able to differentiate them) in order to detect with the cameras, through computer vision techniques, points of interests of it, such as joints and the base. This way, after detecting each marker and saving its coordinates in pixels of every image, triangulation its done so as to get the coordinates in the 3D world, seen from one of the cameras we will tag as the main one. Of course, the cameras were previously calibrated.
With these coordinates, through least squares optimization, we will resolve an equation system that will relate the triangulated coordinates, with the relative coordinates seen from the 'perspective' of our articulated arm. This ecuations are taken from a homogeneus coordinate transformation, from the camere to the base of the arm. So, for every marker, we get 3 equations, that added to the 3+X unknown variables (X degrees of freedom), let us resolve the Yaw, Pitch and Roll relative angles and every degree of freedom (for instance, the degrees a rotative joint is opened) (Xc,Yc,Zc)=F(R,P,Y,Qn).
The resulting images show the value of all of these variables, as well as a graphic representation of how the reference axis of our articulated mechanical system would be seen from the main camera.

THINGS TO KEEP IN MIND:
-Though the place of work had 4 cameras, it was only possible to work with three of them. More cameras means higher precision.
-The OpenCV function that calculates the triangulation, can only work with 1 pair of cameras at a time, so its necessary to average the results.
-These scripts are ready to work with systems of a very particualr configuration, of 2, 3 and 4 degrees of freedom. But it's meant to be a generalized proceedure.
