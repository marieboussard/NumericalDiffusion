function extract_data_stencil2(u, j, sL, sR)
    indices = mod1.(j - sL + 1 : j + sR, size(u, 1))  
    if ndims(u) == 1 
        return view(u,indices)
    else
        return view(u,indices, :)
    end
end

function extract_data_stencil1(u, j, sL, sR)
    Nx = size(u)[1]
    if length(size(u)) > 1
        p = size(u)[2]
    else
        p = 1
    end
    u_short = zeros(eltype(u), (sL+sR,p))
    i = 1
    for k in j-sL+1:j+sR 
        u_short[i,:] = u[mod1(k, Nx),:]
        i+=1
    end
    u_short
end

u = rand(100,3)

j = 4
sL = 8
sR = 8

using Test
@test extract_data_stencil1(u, j, sL, sR) == extract_data_stencil2(u, j, sL, sR)

using ProfileView
@show @allocated extract_data_stencil1(u, j, sL, sR)
@show @allocated extract_data_stencil2(u, j, sL, sR)

@time extract_data_stencil1(u, j, sL, sR)
@time extract_data_stencil2(u, j, sL, sR)