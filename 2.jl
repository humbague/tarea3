function hermitian(n)
    a=rand(n,n)+1im*rand(n,n)
    b=triu(a,1)
    M=b+transpose(conj(b))+Diagonal(rand(n,n))
    return eigvals(M)
end
