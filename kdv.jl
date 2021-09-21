# kdv.jl
##########################
using FFTW, DSP, LinearAlgebra, Plots

# initial data
N=256;
n=128;
M=1000;
u=zeros(ComplexF64, N,M);
v=zeros(ComplexF64, N,M);
z=zeros(ComplexF64, 2*N,M);
w=zeros(ComplexF64, N,M);
U=zeros(ComplexF64, N,N);
V=zeros(ComplexF64, N,M);

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

D1=zeros(ComplexF64, N, N)
D3=zeros(ComplexF64, N, N)
for k=1:n
    D1[k,k]=k-1;
    D1[k+n,k+n]=k-n;
    D3[k,k]=(k-1)^3;
    D3[k+n,k+n]=(k-n)^3;
end

G=im*(8*pi^3*t*N^(-3))*D3;
A=E+G+G^2/2+G^3/6;
B=-3*im*2*pi*t*N^(-1)*(E+G/2+G^2/6+G^3/24)*D1;

v[:,1]=F*u[:,1];
z[:,1]=[v[:,1];v[:,1]];

for k=1:N
    for q=1:N
        w[k,1]=w[k,1]+z[N+k-q+1,1]*z[q,1]/N;
    end
end

for l=2:M
    v[:,l]=A*v[:,l-1]+B*w[:,l-1];
    z[:,l]=[v[:,l];v[:,l]];   
    for k=1:N
        for q=1:N
            w[k,l]=w[k,l]+z[N+k-q+1,l]*z[q,l]/N;
        end
    end
end

U=F'*v/N;

V=zeros(ComplexF64, N,M)
for l=1:M
    V[:,l]=U[:,M-l+1];
end

p=plot(real.(U'), 
   title="KdV flow", 
    st=:heatmap, 
    color=:jet, 
    yaxis="time t", 
    xaxis="position x")
plot!(yticks = ([0 1000;], [0 10]))
plot!(xticks = ([0 256;], [0 256]))

q=surface(real.(V), 
        title="KdV flow", 
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
savefig("kdv.png") 