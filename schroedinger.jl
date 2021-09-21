# schroedinger.jl
##########################
using FFTW, DSP, LinearAlgebra, Plots

# initial data
N=256;
n=128;
M=1000;
u=zeros(ComplexF64, N,M);
for p=65:192
    u[p,1]=exp(2*pi*im*p/128)+exp(2*pi*im*p^2/12800);
end

# setting
E=zeros(ComplexF64, N, N);
for k=1:N
    E[k,k]=1
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

D=zeros(ComplexF64, N, N)
for k=1:n
    D[k,k]=(k-1)^2;
    D[k+n,k+n]=(n-k)^2;
end

G=-im*(4*pi^2*t*N^(-3))*F'*D*F
A=E+G+G^2/2+G^3/6;

for l=2:M
    u[:,l]=A*u[:,l-1];
end

v=zeros(ComplexF64, N,M)
for l=1:M
    v[:,l]=u[:,M-l+1];
end

p1=plot(real.(u'), 
   title="The real part of Schroedinger flow", 
    st=:heatmap, 
    color=:jet, 
    yaxis="time t", 
    xaxis="position x")
plot!(yticks = ([0 1000;], [0 10]))
plot!(xticks = ([0 256;], [0 256]))

p2=plot(imag.(-u'), 
   title="The imaginary part of Schredinger flow", 
    st=:heatmap, 
    color=:jet, 
    yaxis="time t", 
    xaxis="position x")
plot!(yticks = ([0 1000;], [0 10]))
plot!(xticks = ([0 256;], [0 256]))

q1=surface(real.(v), 
        title="The real part of Schroedinger flow", 
        camera=(80,60), 
        color=:jet,
        xaxis="time t", 
        yaxis="position x")
plot!(xticks = ([0 1000;], [10 0]))
plot!(yticks = ([0 256;], [0 256]))

q2=surface(imag.(v), 
        title="The imaginary part of Schroedinger flow", 
        camera=(80,60), 
        color=:jet,
        xaxis="time t", 
        yaxis="position x")
plot!(xticks = ([0 1000;], [10 0]))
plot!(yticks = ([0 256;], [0 256]))
    
plot(p1, q1, p2, q2,  
     layout=(2,2), 
     size=(1600,1300), 
     margin=Plots.Measures.Length(:mm, 10.0))
savefig("schroedinger.png") 