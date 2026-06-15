import math

# configuracion usada en los calculos

# Evita problemas con resultados decimales muy cercanos a cero
TOLERANCIA = 0.0000000001


# funciones para pedir y validar los datos

def pedir_float(mensaje):
    """Pide un numero decimal y valida la entrada"""
    while True:
        try:
            return float(input(mensaje))
        except ValueError:
            print("Entrada no valida. Ingrese un numero.")


def pedir_entero(mensaje, minimo=1):
    """Pide un numero entero mayor o igual al minimo"""
    while True:
        try:
            numero = int(input(mensaje))
            if numero >= minimo:
                return numero
            print("Debe ingresar un valor mayor o igual que {}.".format(minimo))
        except ValueError:
            print("Entrada no valida. Ingrese un numero entero.")


# funciones para mostrar tablas matrices y resultados

def imprimir_matriz(matriz):
    """Muestra una matriz aumentada con cuatro decimales"""
    # La ultima columna contiene los resultados del sistema
    for fila in matriz:
        texto = "[ "
        for i in range(len(fila)):
            valor = 0.0 if abs(fila[i]) < TOLERANCIA else fila[i]
            if i == len(fila) - 1:
                texto += "| "
            texto += "{:10.4f} ".format(valor)
        print(texto + "]")
    print()


def imprimir_tabla(titulos, filas, nombre="TABLA DE MINIMOS CUADRADOS"):
    """Muestra los datos organizados en una tabla"""
    print("\n" + nombre)
    for titulo in titulos:
        print("{:>14}".format(titulo), end="")
    print()
    for fila in filas:
        for valor in fila:
            print("{:14.4f}".format(valor), end="")
        print()


def mostrar_sumatorias(nombres, valores, cantidad):
    """Muestra las sumatorias y la cantidad de datos"""
    print("\nSUMATORIAS")
    for i in range(len(nombres)):
        print("Suma {:<6} = {:.4f}".format(nombres[i], valores[i]))
    print("n           = {}".format(cantidad))


def mostrar_sistema(matriz):
    """Muestra las ecuaciones normales del sistema"""
    # Cada fila representa una ecuacion con a0 a1 y a2
    print("\nSISTEMA DE ECUACIONES NORMALES")
    cantidad = len(matriz)
    for fila in matriz:
        ecuacion = ""
        for i in range(cantidad):
            if i > 0:
                ecuacion += " + "
            ecuacion += "{:.4f}a{}".format(fila[i], i)
        print(ecuacion + " = {:.4f}".format(fila[-1]))


# funciones para resolver los sistemas de ecuaciones

def buscar_pivote(matriz, pivote):
    """Busca otra fila cuando el pivote es cero"""
    # Se cambia la fila para poder continuar con el procedimiento
    if abs(matriz[pivote][pivote]) >= TOLERANCIA:
        return True
    for fila in range(pivote + 1, len(matriz)):
        if abs(matriz[fila][pivote]) >= TOLERANCIA:
            matriz[pivote], matriz[fila] = matriz[fila], matriz[pivote]
            print("Se intercambia F{} con F{}.".format(pivote + 1, fila + 1))
            imprimir_matriz(matriz)
            return True
    return False


def eliminacion_gaussiana(matriz):
    """Resuelve el sistema por eliminacion gaussiana"""
    # En este metodo se hacen ceros debajo de cada pivote
    # Al final se usa sustitucion hacia atras para encontrar a0 a1 y a2

    # Se trabaja con una copia para conservar la matriz original
    m = [fila[:] for fila in matriz]
    cantidad = len(m)
    print("\nMATRIZ AUMENTADA INICIAL")
    imprimir_matriz(m)

    for pivote in range(cantidad):
        # Se comprueba que el pivote no sea cero
        if not buscar_pivote(m, pivote):
            print("El sistema no tiene una solucion unica.")
            return None

        # Se divide la fila para convertir el pivote en 1
        divisor = m[pivote][pivote]
        for columna in range(pivote, cantidad + 1):
            m[pivote][columna] /= divisor
        print("Convertir el pivote {} en 1.".format(pivote + 1))
        imprimir_matriz(m)

        # Se resta un multiplo de la fila pivote a las filas de abajo
        # De esta forma los numeros debajo del pivote quedan en cero
        for fila in range(pivote + 1, cantidad):
            factor = m[fila][pivote]
            for columna in range(pivote, cantidad + 1):
                m[fila][columna] -= factor * m[pivote][columna]
        if pivote < cantidad - 1:
            print("Hacer ceros debajo del pivote {}.".format(pivote + 1))
            imprimir_matriz(m)

    # Se empieza por la ultima ecuacion porque tiene menos incognitas
    # Luego se sube fila por fila hasta encontrar todos los coeficientes
    soluciones = [0.0] * cantidad
    print("SUSTITUCION HACIA ATRAS")
    for fila in range(cantidad - 1, -1, -1):
        suma = 0.0
        for columna in range(fila + 1, cantidad):
            suma += m[fila][columna] * soluciones[columna]
        soluciones[fila] = m[fila][-1] - suma
        print("a{} = {:.4f}".format(fila, soluciones[fila]))
    return soluciones


def gauss_jordan(matriz):
    """Resuelve el sistema por Gauss Jordan"""
    # Este metodo busca convertir la parte izquierda en una matriz identidad

    # Se trabaja con una copia de la matriz
    m = [fila[:] for fila in matriz]
    cantidad = len(m)
    print("\nMATRIZ AUMENTADA INICIAL")
    imprimir_matriz(m)

    for pivote in range(cantidad):
        if not buscar_pivote(m, pivote):
            print("El sistema no tiene una solucion unica.")
            return None

        # Se convierte el pivote en 1
        divisor = m[pivote][pivote]
        for columna in range(cantidad + 1):
            m[pivote][columna] /= divisor
        print("Paso {}: convertir el pivote en 1.".format(pivote + 1))
        imprimir_matriz(m)

        # Se resta la fila pivote a todas las otras filas
        # Asi se forman ceros arriba y abajo del pivote
        for fila in range(cantidad):
            if fila != pivote:
                factor = m[fila][pivote]
                for columna in range(cantidad + 1):
                    m[fila][columna] -= factor * m[pivote][columna]
        print("Hacer ceros arriba y abajo del pivote.")
        imprimir_matriz(m)

    print("MATRIZ FINAL")
    imprimir_matriz(m)
    return [m[i][-1] for i in range(cantidad)]


# funciones para calcular y mostrar la calidad del ajuste

def calcular_r2(y_reales, y_estimados):
    """Calcula el coeficiente de determinacion r2"""
    # Se calcula el promedio de los valores reales
    promedio = sum(y_reales) / len(y_reales)
    error = 0.0
    total = 0.0
    for i in range(len(y_reales)):
        # El error compara el valor real con el valor estimado
        # El total compara el valor real con el promedio
        error += (y_reales[i] - y_estimados[i]) ** 2
        total += (y_reales[i] - promedio) ** 2

    # No se puede dividir si la variacion total es cero
    if abs(total) < TOLERANCIA:
        return None, promedio, error, total

    # Se aplica la formula vista en clase
    # r2 = 1 - suma de errores entre suma de variacion total
    r2 = 1 - error / total
    if abs(r2) < TOLERANCIA:
        r2 = 0.0
    elif abs(r2 - 1) < TOLERANCIA:
        r2 = 1.0
    return r2, promedio, error, total


def interpretar_r2(r2):
    """Indica la calidad del ajuste obtenido"""
    if r2 < 0.20:
        return "Muy debil"
    elif r2 < 0.40:
        return "Debil"
    elif r2 < 0.60:
        return "Moderada"
    elif r2 < 0.80:
        return "Fuerte"
    return "Muy fuerte"


def mostrar_r2(y_reales, y_estimados, pendiente=None):
    """Muestra los calculos y la interpretacion de r2"""
    r2, promedio, error, total = calcular_r2(y_reales, y_estimados)

    # La tabla permite ver el error de cada dato por separado
    imprimir_tabla(
        ["y", "y estimada", "(y-y est.)^2", "(y-prom.)^2"],
        [
            [
                y_reales[i],
                y_estimados[i],
                (y_reales[i] - y_estimados[i]) ** 2,
                (y_reales[i] - promedio) ** 2,
            ]
            for i in range(len(y_reales))
        ],
        "TABLA PARA r^2",
    )
    print("\ny promedio = {:.4f}".format(promedio))
    print("Suma (y - y estimada)^2 = {:.4f}".format(error))
    print("Suma (y - y promedio)^2 = {:.4f}".format(total))

    if r2 is None:
        print("No se puede calcular r^2 porque todos los valores de y son iguales.")
        return

    print("\nCOEFICIENTE DE DETERMINACION")
    print("r^2 = 1 - ({:.4f} / {:.4f}) = {:.4f}".format(error, total, r2))
    if pendiente is not None:
        # El signo de r depende del signo de la pendiente
        r = math.sqrt(max(0, r2))
        if pendiente < 0:
            r = -r
        print("\nCOEFICIENTE DE CORRELACION")
        print("r = {:.4f}".format(r))

    print("\nINTERPRETACION")
    print("El modelo explica aproximadamente el {:.4f}% de la variabilidad de y.".format(r2 * 100))
    print("Interpretacion del ajuste: {}.".format(interpretar_r2(r2)))


# funciones para guardar datos y formar la ecuacion

def obtener_datos(nombres, minimo):
    """Pide los datos necesarios para realizar la regresion"""
    cantidad = pedir_entero("Ingrese la cantidad de datos: ", minimo)
    # Se crea una lista para cada variable
    columnas = [[] for nombre in nombres]

    # Se guardan los datos una fila a la vez
    for i in range(cantidad):
        print("\nDato {}".format(i + 1))
        for j in range(len(nombres)):
            columnas[j].append(pedir_float("Ingrese {}: ".format(nombres[j])))
    return columnas


def mostrar_resultado(coeficientes, variables):
    """Muestra los coeficientes y la ecuacion final"""
    print("\nCOEFICIENTES")
    for i in range(len(coeficientes)):
        print("a{} = {:.4f}".format(i, coeficientes[i]))

    # Se coloca el signo correcto de cada coeficiente
    texto = "y estimada = {:.4f}".format(coeficientes[0])
    for i in range(1, len(coeficientes)):
        signo = "+" if coeficientes[i] >= 0 else "-"
        texto += " {} {:.4f}{}".format(signo, abs(coeficientes[i]), variables[i - 1])
    print("\nECUACION DEL MODELO")
    print(texto)


# metodos de regresion

def regresion_lineal():
    """Realiza la regresion lineal simple"""
    print("\n" + "-" * 60)
    print("REGRESION LINEAL SIMPLE")
    print("Modelo: y estimada = a0 + a1x")
    print("-" * 60)

    x, y = obtener_datos(["x", "y"], 2)

    # Se calculan las columnas x y x al cuadrado y xy
    filas = [[x[i], y[i], x[i] ** 2, x[i] * y[i]] for i in range(len(x))]
    imprimir_tabla(["x", "y", "x^2", "xy"], filas)

    # Se suman las columnas de la tabla
    sx, sy, sx2, sxy = [sum(fila[i] for fila in filas) for i in range(4)]
    mostrar_sumatorias(["x", "y", "x^2", "xy"], [sx, sy, sx2, sxy], len(x))

    # Estas son las dos ecuaciones normales del metodo
    # n*a0 + Sx*a1 = Sy
    # Sx*a0 + Sx2*a1 = Sxy
    matriz = [[len(x), sx, sy], [sx, sx2, sxy]]
    mostrar_sistema(matriz)

    # El determinante permite saber si el sistema tiene una sola solucion
    # Si es cero no se pueden calcular coeficientes unicos
    determinante = len(x) * sx2 - sx ** 2
    if abs(determinante) < TOLERANCIA:
        print("No existe una solucion unica.")
        return

    a0 = (sy * sx2 - sx * sxy) / determinante
    a1 = (len(x) * sxy - sx * sy) / determinante
    print("\nDeterminante = {:.4f}".format(determinante))
    mostrar_resultado([a0, a1], ["x"])

    # Se usa el modelo para realizar el pronostico
    xp = pedir_float("\nIngrese x para pronosticar: ")
    print("\nPRONOSTICO")
    print("Para x = {:.4f}, y estimada = {:.4f}".format(xp, a0 + a1 * xp))
    # Se calculan los valores estimados para obtener r2
    mostrar_r2(y, [a0 + a1 * valor for valor in x], a1)


def regresion_multiple():
    """Realiza la regresion lineal multiple"""
    print("\n" + "-" * 60)
    print("REGRESION LINEAL MULTIPLE")
    print("Modelo: y estimada = a0 + a1x1 + a2x2")
    print("-" * 60)

    x1, x2, y = obtener_datos(["x1", "x2", "y"], 3)
    # Se forman las columnas necesarias para los calculos
    filas = []
    for i in range(len(y)):
        filas.append([
            x1[i], x2[i], y[i], x1[i] * x2[i], x1[i] ** 2,
            x2[i] ** 2, x1[i] * y[i], x2[i] * y[i],
        ])
    nombres = ["x1", "x2", "y", "x1x2", "x1^2", "x2^2", "x1y", "x2y"]
    imprimir_tabla(nombres, filas)

    # Se suman todas las columnas
    sumas = [sum(fila[i] for fila in filas) for i in range(8)]
    mostrar_sumatorias(nombres, sumas, len(y))
    sx1, sx2, sy, sx1x2, sx1_2, sx2_2, sx1y, sx2y = sumas

    # Cada fila de la matriz representa una ecuacion normal
    # La ultima columna contiene los resultados de las ecuaciones
    matriz = [
        [len(y), sx1, sx2, sy],
        [sx1, sx1_2, sx1x2, sx1y],
        [sx2, sx1x2, sx2_2, sx2y],
    ]
    mostrar_sistema(matriz)
    # Se resuelve la matriz por eliminacion gaussiana
    coeficientes = eliminacion_gaussiana(matriz)
    if coeficientes is None:
        return
    mostrar_resultado(coeficientes, ["x1", "x2"])

    xp1 = pedir_float("\nIngrese x1 para pronosticar: ")
    xp2 = pedir_float("Ingrese x2 para pronosticar: ")
    # Se sustituyen x1 y x2 para hacer el pronostico
    pronostico = coeficientes[0] + coeficientes[1] * xp1 + coeficientes[2] * xp2
    print("\nPRONOSTICO")
    print("Para x1 = {:.4f} y x2 = {:.4f}, y estimada = {:.4f}".format(xp1, xp2, pronostico))

    # Se calculan los valores estimados para obtener r2
    estimados = [
        coeficientes[0] + coeficientes[1] * x1[i] + coeficientes[2] * x2[i]
        for i in range(len(y))
    ]
    mostrar_r2(y, estimados)


def regresion_polinomial():
    """Realiza la regresion polinomial de segundo grado"""
    print("\n" + "-" * 60)
    print("REGRESION POLINOMIAL DE 2DO GRADO")
    print("Modelo: y estimada = a0 + a1x + a2x^2")
    print("-" * 60)

    x, y = obtener_datos(["x", "y"], 3)
    # Se necesitan potencias hasta x4 para formar las ecuaciones normales
    # Tambien se calculan xy y x2y
    filas = []
    for i in range(len(x)):
        filas.append([
            x[i], y[i], x[i] ** 2, x[i] ** 3, x[i] ** 4,
            x[i] * y[i], (x[i] ** 2) * y[i],
        ])
    nombres = ["x", "y", "x^2", "x^3", "x^4", "xy", "x^2y"]
    imprimir_tabla(nombres, filas)

    # Se suman las columnas de la tabla
    sumas = [sum(fila[i] for fila in filas) for i in range(7)]
    mostrar_sumatorias(nombres, sumas, len(x))
    sx, sy, sx2, sx3, sx4, sxy, sx2y = sumas

    # Cada fila de la matriz representa una ecuacion normal
    # La ultima columna contiene los resultados de las ecuaciones
    matriz = [
        [len(x), sx, sx2, sy],
        [sx, sx2, sx3, sxy],
        [sx2, sx3, sx4, sx2y],
    ]
    mostrar_sistema(matriz)
    # Se resuelve la matriz por Gauss Jordan
    coeficientes = gauss_jordan(matriz)
    if coeficientes is None:
        return
    mostrar_resultado(coeficientes, ["x", "x^2"])

    xp = pedir_float("\nIngrese x para pronosticar: ")
    # Se sustituye x para hacer el pronostico
    pronostico = coeficientes[0] + coeficientes[1] * xp + coeficientes[2] * xp ** 2
    print("\nPRONOSTICO")
    print("Para x = {:.4f}, y estimada = {:.4f}".format(xp, pronostico))

    # Se calculan los valores estimados para obtener r2
    estimados = [
        coeficientes[0] + coeficientes[1] * valor + coeficientes[2] * valor ** 2
        for valor in x
    ]
    mostrar_r2(y, estimados)


# menu principal

def menu_principal():
    """Muestra el menu principal del programa"""
    # El menu se repite hasta seleccionar la opcion 4
    while True:
        print("\n" + "=" * 40)
        print("        PROYECTO DE REGRESION")
        print("=" * 40)
        print("1. Regresion Lineal Simple")
        print("2. Regresion Lineal Multiple")
        print("3. Regresion Polinomial de 2do Grado")
        print("4. Salir")
        opcion = input("Seleccione una opcion: ").strip()

        if opcion == "1":
            regresion_lineal()
        elif opcion == "2":
            regresion_multiple()
        elif opcion == "3":
            regresion_polinomial()
        elif opcion == "4":
            print("\nPrograma finalizado.")
            break
        else:
            print("Opcion no valida.")


# Inicia el programa
if __name__ == "__main__":
    menu_principal()
