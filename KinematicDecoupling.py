from sympy import*
from sympy.matrices import*
init_printing(use_latex = true)
theta1,theta2,theta3,theta4,theta5,theta6 = symbols("theta_1, theta_2, theta_3, theta_4, theta_5, theta_6")

#// Funcion para obtener matriz Denavit-Hartenberg
def MDH(ai, alphai, di, thetai):
    thetaic = thetai
    alphaic = alphai

    M = Matrix([
        [cos(thetaic), -sin(thetaic)*cos(alphaic),  sin(thetaic)*sin(alphaic), ai*cos(thetaic)],
        [sin(thetaic),  cos(thetaic)*cos(alphaic), -cos(thetaic)*sin(alphaic), ai*sin(thetaic)],
        [0,             sin(alphaic),               cos(alphaic),              di],
        [0,             0,                          0,                         1]
    ])
    return M

#// Funcion desacoplo cinematico
def KinematicDecoupling(H):
    """
    Calcular la cinematica inversa del robot ABB IRB 140 por medio del metodo 
    desacoplo cinematico
    
    Parametros
    --------------------
    H : Matriz de transforcion homogenea que describe la posicion y orientacion 
    deseads del extremo del robot ABB IRB 140
    
    Devoluciones
    --------------------
    Angulos: La funcion devuelve los angulos requeridos de cada una de las
    articulaciones del robot angulos = {θ1, θ2, θ3, θ4, θ5, θ6} en radianes
    
    HOUT : La funcion devuelve la matriz evaluada con los angulos de las 
    articulaciones obtenidas
    
    Ejemplos:
    --------------------
    Matriz de transformacion homogenea que recibe la funcion :
    H = Matrix([
    [-sqrt(3)/4, -1/4, sqrt(3)/2, 500],
    [1/2, -sqrt(3)/2, 0, 100],
    [3/4, sqrt(3)/4, 1/2, 1000],
    [0, 0, 0, 1]
    ])
    
    [HOUT, angles] = KinematicDecoupling(H)
    """
    #// Parametros de DH para el robot ABB IRB 140
    T01 = MDH(70, pi/2, 352, theta1)
    T12 = MDH(360, 0, 0, theta2)
    T23 = MDH(0, pi/2, 0, theta3)
    T34 = MDH(0, -pi/2, 380, theta4)
    T45 = MDH(0, pi/2, 0, theta5)
    T56 = MDH(0, 0, 65, theta6)
    #// Matriz que define el sistema 6 en 0
    T06 = T01*T12*T23*T34*T45*T56
    #// Obtener matriz de rotacion de 3 en 0
    T03 = T01*T12*T23
    R03 = simplify(T03[:3,:3])
    #// Obtener matriz de rotacion de 6 en 0
    R06  = simplify(H[:3,:3])
    
    #// Metodo geometrico para obtener los angulos de las 3 primeras articulaciones
    [Xm, Ym, Zm] = H[:3,3].T - 65*H[:3,2].T

    ta1 = (atan(Ym/Xm))

    a = sqrt(Xm**2 + Ym**2)
    b = Zm - 352
    c = a - 70
    r = sqrt(b**2 + c**2)

    betha = (atan(b/c))
    k = ((360**2 + r**2 - 380**2)/(2*r*360))
    gamma = (atan2(sqrt(1 - k**2),k))
    ta2 = (betha - gamma)

    k1 = ((r**2 - 360**2 - 380**2)/(2*360*380))
    alpha = (atan2(sqrt(1 - k1**2),k1))
    ta3 = (alpha + pi/2)
    
    t1 = ta1
    t2 = ta2
    t3 = ta3
    #// Sustituir valores obtenidos de las 3 primeras articulaciones en la matriz de rotacion 6 en 3
    NR03 = R03.subs(theta1,t1).subs(theta2,t2).subs(theta3,t3)
    NR03T = NR03.T
    R36 = (NR03T)*(R06)
    R36.evalf()
    #// Matriz de rotacion 6 en 3 obtenida de la cinematica directa
    T36 = T34*T45*T56
    #// Obtner los valores de las articulaciones que proporcionan la orientacion igualando matrices y despejando
    t4 = (atan2(1*R36[1,2], 1*R36[0,2]))
    t5 = (atan2(sqrt(1 - R36[2,2]**2), R36[2,2]))
    t6 = (atan2(1*R36[2,1], -R36[2,0]))
    
    NT01 = MDH(70, pi/2, 352, t1.evalf())
    NT12 = MDH(360, 0, 0, t2.evalf())
    NT23 = MDH(0, pi/2, 0, t3.evalf())
    NT34 = MDH(0, -pi/2, 380, t4.evalf())
    NT45 = MDH(0, pi/2, 0, t5.evalf())
    NT56 = MDH(0, 0, 65, t6.evalf())
    
    HOUT = NT01*NT12*NT23*NT34*NT45*NT56
    Angulos = [t1.evalf(), t2.evalf(), t3.evalf(), t4.evalf(), t5.evalf(), t6.evalf()]
    
    return HOUT, Angulos
