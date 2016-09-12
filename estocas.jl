using Base.Test

function stoca(n,m)
    A=zeros(n,m)
    for i in 1:n
        v=rand(m)
        A[i,:]=v/sum(v)
    end
    return A
end

function prueba(A)
    for i in 1:size(A,1)
        if abs(sum(A,2)[i]-1.0)>1e-5
            return false
        end
    end
    return true
end

@test prueba(stoca(3,2))
