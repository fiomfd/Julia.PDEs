# lkdv.jl
##########################
using FFTW, DSP, LinearAlgebra, Plots

# initial data
N=256;
n=128;
M=1000;
u=zeros(ComplexF64, N,M);
v=zeros(ComplexF64, N,M);
for p=33:96
    u[p,1]=1/4;
    u[p+128,1]=1/8;
end

# setting
E=zeros(ComplexF64, N, N);
for k=1:N
    E[k,k]=1;
end

t=1/100;

# sqrt(M) times Fourier matrix
omega=exp(im*2*pi/N);
F=zeros(ComplexF64, N, N);
for k=1:N
    for l=1:N
    F[k,l]=omega^(-(k-1)*(l-1));
    end
end

v[:,1]=F*u[:,1];
D1=zeros(ComplexF64, N, N)
D3=zeros(ComplexF64, N, N)
for k=1:n
    D1[k,k]=k-1;
    D1[k+n,k+n]=k-n;
    D3[k,k]=(k-1)^3;
    D3[k+n,k+n]=(k-n)^3;
end

G=F'*(im*8*pi^3*t*N^(-4))*D3*F;
A=E+G+G^2/2+G^3/6;

for l=2:M
    u[:,l]=A*u[:,l-1];
end

v=zeros(ComplexF64, N,M)
for l=1:M
    v[:,l]=u[:,M-l+1];
end

p=plot(real.(u'), 
       title="Linearlized KdV flow", 
       st=:heatmap, 
       color=:jet, 
       yaxis="time t", 
       xaxis="position x")
plot!(yticks = ([0 1000;], [0 10]))
plot!(xticks = ([0 256;], [0 256]))

q=surface(real.(v), 
        title="Linearlized KdV flow", 
        camera=(80,60), 
        color=:jet,
        xaxis="time t", 
        yaxis="position x")
plot!(xticks = ([0 1000;], [10 0]))
plot!(yticks = ([0 256;], [0 256]))
    
plot(p, q,  
     layout=(1,2), 
     size=(1600,600), 
     margin=Plots.Measures.Length(:mm, 10.0))
savefig("lkdv.png") 