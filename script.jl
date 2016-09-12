using Base.Test

function hermitian(n::Int)
    a=rand(n,n)+1im*rand(n,n)
    b=triu(a,1)
    M=b+transpose(conj(b))+Diagonal(rand(n,n))
    return M
end

function daga(state::Array{Complex{Float64},1})
    return transpose(conj(state))
end

function proyeccion(dim::Int)
    M=hermitian(dim)
    A=zeros(dim,dim)
    vecs=eigvecs(M)
    for i in 1:dim
        A+=kron(vecs[:,i],daga(vecs[:,i]))
    end
    # chop(A)
    return A
end

#La función prueba probará si la matriz obtenida a través de la función proyeccion es unitaria
function prueba(dim::Int)
    MM=hermitian(dim)
    M=proyeccion(dim)
    A=eye(dim)-M
    b=0
    for i in A
        b+=abs(i)
    end
    if b>1e-5
        return false
    end
    return true
end
    
@test ishermitian(hermitian(3))    
@test prueba(3)
