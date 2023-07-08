import numpy as np

import math

"""
O pacote numpy é utilizado para resolver o sistema linear obtido, já o pacote
math é utilizado para definir a função que corresponde à solução analítica
"""

L=math.pi
T=20
c=1
b=0.2
m=350

"""
L é o comprimento da corda
T é o tempo máximo
c,b são as constantes da EDP
m é o índice do somatório da solução analíca
"""

def f(x):
    return 1-(2*abs(x-0.5*L))/L

def sin(x):
    return(math.sin(x))

def cos(x):
    return(math.cos(x))

def exp(x):
    return(math.exp(x))

def sqrt(x):
    return(math.sqrt(x))

def g(x,t,n):
    return (8/(n*n*math.pi*math.pi)*sin(0.5*n*math.pi)*sin(n*math.pi*x/L)*exp(-b*t)
            *(cos(t*sqrt(n*n*math.pi*math.pi/(L*L)*c*c-b*b))+(b*L)
            /(sqrt(n*n*math.pi*math.pi*c*c-L*L*b*b))
            *sin(t*sqrt(n*n*math.pi*math.pi/(L*L)*c*c-b*b))))

def u(x,t):
    y=0
    for n in range(1,m):
       y=y+g(x,t,n)
    return y

"""
f(x) é a função descrita nas condições de contorno
g(x,t,n) é a função que associa para cada valor de n a espressão no somatório
da solução analítica
u(x,t) soma as expressões de g(x,t,n) e é, portanto, a solução analítica de fato 
"""

h=L/6
k=T/40

l=c*k/h

r=float(L/h-1)
s=float(T/k)

"""
h é a distância horizontal entre os pontos da malha
k é a distância vertical entre os pontos da malha
l é a constante lambda associada com a convergência do método
r é a quantidade de pontos internos em cada linha da malha
s é a quantidade de linhas na malha
"""

v=[]

"""
O seguinte loop cria a matriz dos coeficientes do problema.
Os pontos estão organizados por linha, de modo que todos os pontos na primeira
posição de cada linha (1, r+1, 2r+1, ...) não tem vizinhos à esquerda e os 
pontos na última posição de cada linha (r, 2r, 3r, ...) não tem vizinhos à 
direita. 
O código então associada a cada ponto u_{ij} a expressão l^2u_{i-1,j}+
2(1-l^2)u_{ij}+l^2u_{i+1,j}+(bk-1)u_{i,j-1}-(bk+1)u_{i,j+1}.
As últimas r entradas correspondem às equações obtidas para os pontos da linha
1 a partir dos pontos da linha 0 e, portanto, tem apenas um elemento em cada
coluna.
"""

for i in range(1,int(r*s+1)):
    vi=[]
    if i<(int(r*(s-1)+1)):
        if (i%r)==0:
            for j in range(1,int(r*s+1)):
                if j==(i-r):
                    vi.append(b*k-1)
                elif j==(i-1):
                    vi.append(l*l)
                elif j==i:
                    vi.append(2*(1-l*l))
                elif j==(i+r):
                    vi.append(-1-b*k)
                else:
                    vi.append(0)
        elif (i%r)==1:
            for j in range(1,int(r*s+1)):
                if j==(i-r):
                    vi.append(b*k-1)
                elif j==i:
                    vi.append(2*(1-l*l))
                elif j==(i+1):
                    vi.append(l*l)
                elif j==(i+r):
                    vi.append(-1-b*k)
                else:
                    vi.append(0)
        else:
            for j in range(1,int(r*s+1)):
                if j==(i-r):
                    vi.append(b*k-1)
                elif j==(i-1):
                    vi.append(l*l)
                elif j==i:
                    vi.append(2*(1-l*l))
                elif j==(i+1):
                    vi.append(l*l)
                elif j==(i+r):
                    vi.append(-1-b*k)
                else:
                    vi.append(0)
    else:
        for j in range(1,int(r*s+1)):
            if j<(r+1):
                if ((j-i)%r)==0:
                    vi.append(1)
                else:
                    vi.append(0)
            else:
                vi.append(0)
    v.append(vi)    

A=np.array(v)

"""
O vetor B é o vetor dos termos independentes
"""

B=[]

for i in range(1,int(r*s+1)):
    if i<=r:
        B.append((1-b*k)*f(i*h))
    elif i<=(r*(s-1)):
        B.append(0)
    elif i<=(r*s-1):
        B.append(0.5*l*l*(f(float(((i-1)%r)*h))+f(float((i%r+1))*h))
        +(1-l*l)*(f(float((i%r)*h))))
    else:
        B.append(0.5*l*l*(f((float(r)-1)*h))+(1-l*l)*(f(r*h)))

X = np.linalg.inv(A).dot(B)

W = np.array_split(X,int(s))

"""
O vetor X é a solução do sistema.
W é o vetor X com as entradas separadas por linha.
"""

Y=[]  

for i in range(1,int(r+1)):
    for j in range(1,int(s+1)):
        Y.append(u(i*h,j*k))

Y1=np.array_split(Y,r)

"""
O vetor Y contém o valor da solução analítica em cada ponto. Os valores de Y
estão organizados por coluna, portanto, o vetor Y1 separa os valores em colunas e 
o vetor Z é criado para reorganizá-lo por linhas da mesma forma que o vetor W,. 
"""

Z=[]

for i in range(0,len(Y1[1])):
    vi=[]
    for j in range(0,len(Y1)):
        vi.append(Y1[j][i])
    Z.append(vi)

"""
O loop a seguir cria o código em latex que produz a tabela com os valores 
aproximados de alguns dos pontos da malha.
"""

texto = "\\begin{table}[H]\\centering\\renewcommand{\\arraystretch}{1.5}\\begin{tabular}{"

for i in range(0,int(r)+1):
    texto=texto+"c "

texto=texto+"}\hline Tempo"

for i in range(1,int(r)+1):
    texto=texto+" & x="+str(round(i*h,3))

texto=texto+" \\\ \hline "

for i in range(0,len(W)):
    if i%2==1:
        tex = ""
        for j in range(0,len(W[i])):
            if j==0:
                tex=tex+str(round((i+1)*k,3))+" & "+str(round(W[i][j],4))+" & "
            elif j!=(len(W[i])-1):
                tex=tex+str(round(W[i][j],4))+" & "
            else:
                tex=tex+str(round(W[i][j],4))
        texto = texto+tex+" \\\ "
    
texto=texto+"\\end{tabular}\\end{table}"

"""
O vetor E contém os erros associados a cada valor presente na tabela obtida
com o loop anterior.
"""

E=[]

for i in range(0,len(W)):
    vi=[]
    for j in range(0,len(W[i])):
        vi.append((W[i][j]-Z[i][j])/W[i][j])
    E.append(vi)
    
"""
O loop a seguir cria o código em latex que produz a tabela com os erros
associados aos valores aproximados da tabela anterior.
"""

texto1 = "\\begin{table}[H]\\centering\\renewcommand{\\arraystretch}{1.5}\\begin{tabular}{"

for i in range(0,int(r)+1):
    texto1=texto1+"c "

texto1=texto1+"}\hline Tempo"

for i in range(1,int(r)+1):
    texto1=texto1+" & x="+str(round(i*h,3))

texto1=texto1+" \\\ \hline "

for i in range(0,len(E)):
    if i%2==1:
        tex = ""
        for j in range(0,len(E[i])):
            if j==0:
                tex=tex+str(round((i+1)*k,3))+" & "+str(round(E[i][j],4))+" & "
            elif j!=(len(E[i])-1):
                tex=tex+str(round(E[i][j],4))+" & "
            else:
                tex=tex+str(round(E[i][j],4))
        texto1 = texto1+tex+" \\\ "
    
texto1=texto1+"\\end{tabular}\\end{table}"
    
print("Valores aproximados:")
print("")
print(texto)
print("")
print("Erro relativo:")
print("")
print(texto1)