import numpy as np
import matplotlib.pyplot as plt
from math import pi


def main():

    y = np.array([30.8, 31.0, 30.5, 30.0, 29.3, 27.7, 26.4, 23.8, 21.3, 18.7, 17.2, 15.8, 12.7, 10.5, 9.9, 8.7, 7.7, 6.7, 6.2, 6.4, 6.6, 7.1, 8.1, 9.5, 10.3, 12.2, 13.9, 15.1, 18.0, 19.8, 21.2, 22.4, 25.1, 26.0, 28.0, 28.9, 29.6])
    aprox_0 = 10 ** (-5)
    erros_ang = np.arange(-0.5, 0.5, 0.001)
    #erros_ang = np.arange(-0.5, 0.5, 0.00001)
    erros_ang = np.array(erros_ang)

    #
    desvios = []

    for erro_ang in erros_ang:
        (cos,ang) = cos_ang(erro_ang,aprox_0)
        (m, n) = achar_regressao(cos, y, 37)
        desvio = desvio_padrao(cos, y, 37, m, n)
        desvios.append(desvio)
        if desvio<=min(desvios):
            erro_ang_min=erro_ang
            (m_min,n_min)=(m,n)
        #erro_ang=-0.26888

    desvio_min = min(desvios)
    plt.xlabel(r'$\theta$', fontsize=14)
    plt.ylabel('Desvio-Padrão de Resíduos', fontsize=14)
    plt.text(erro_ang_min-0.15,desvio_min-0.1, 'Erro angular de maior precisão',fontsize=8)
    plt.plot(erros_ang, desvios)
    plt.grid()
    plt.savefig('DesvioXErro.png')
    plt.show()
    #

    print(desvio_min,erro_ang_min,"desvio=",desvio_min)


    (cos_min, ang_min) = cos_ang(erro_ang_min, aprox_0)
    plt.scatter(cos_min,y)
    plt.plot([0,1],[n_min,m_min+n_min])
    titulo1='Intensidade X cos^2 com '+str(round(erro_ang_min,7))+'° de correção'
    plt.title(titulo1)
    text_erro='Desvio-Padrão de Resíduos: '+str(round(desvio_min,3))
    plt.text(0.3,11/2,text_erro)
    text_i='I(0) = '+str(round(n_min,3))
    plt.text(0.01, n_min-1.4, text_i)
    plt.axis([0, 1, 0, 33])
    plt.xlabel('cos^2'r'$(\theta$)', fontsize=14)
    plt.ylabel('Intensidade', fontsize=14)
    plt.grid()
    plt.savefig('Cos_erro.png')
    plt.show()


    (cos0, ang0) = cos_ang(0, aprox_0)
    (m0, n0) = achar_regressao(cos0, y, 37)
    desvio0 = desvio_padrao(cos0, y, 37, m0, n0)
    plt.scatter(cos0, y)
    plt.plot([0, 1], [n0, m0 + n0])
    titulo2 = 'Intensidade X cos^2 sem correção'
    plt.title(titulo2)
    text_erro2 = 'Desvio-Padrão de Resíduos: ' + str(round(desvio0, 3))
    plt.text(0.3, 11 / 2, text_erro2)
    plt.axis([0, 1, 0, 33])
    plt.xlabel('cos^2'r'$(\theta$)', fontsize=14)
    plt.ylabel('Intensidade', fontsize=14)
    plt.grid()
    plt.savefig('Cos_sem.png')
    plt.show()



    a = np.linspace(0, pi)
    plt.plot(a*180/pi, n_min+np.cos(a)*np.cos(a)*m_min)
    plt.plot(180 * ang_min / pi, y,linestyle='solid')
    plt.axis([0, 170, 0, 33])
    titulo3='Intensidade X ângulo corrigido de '+str(round(erro_ang_min,7))+'°'
    plt.title(titulo3)
    plt.xlabel(r'$\theta$', fontsize=14)
    plt.ylabel('Intensidade', fontsize=14)
    plt.grid(True)
    plt.savefig('cos_ang.png')
    plt.show()

    #poderia ter usado uma função de plotagem

def achar_regressao(x, y, n):
    # y=mx+n
    prod_xy = 0
    soma_x = 0
    soma_x2 = 0
    soma_y2 = 0
    soma_y = 0
    numero = 0
    for i in range(n):
        prod_xy += x[i] * y[i]
        soma_x += x[i]
        soma_y += y[i]
        soma_x2 += x[i] * x[i]
        soma_y2 += y[i] * y[i]
        numero += 1
    m = (numero * prod_xy - soma_x * soma_y) / (numero * soma_x2 - soma_x * soma_x)
    n = (soma_y - m * soma_x) / numero
    return (m, n)

def cos_ang(erro,aprox_0):
    cos = np.zeros(37)
    for i in range(37):
        cos[i] = (i * (5 + erro) * np.pi / 180)
    ang=cos
    cos=np.cos(cos)*np.cos(cos)
    for j in range(36):
        if (cos[j] < aprox_0):
            cos[j] = 0
    return (cos,ang)

def desvio_padrao(cos, y, n, m_c, n_c):
    erro = 0
    for i in range(n):
        erro += (y[i] - (m_c) * cos[i] - n_c) * (y[i] - (m_c) * cos[i] - n_c)
    return np.sqrt(erro / (n - 2))

main()