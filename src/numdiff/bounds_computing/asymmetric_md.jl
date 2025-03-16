function utilde!(::SymmetricMD, estimator::Estimator, j::Int)# compute ̃uᵢʲ⁺¹/²
    @unpack uinit = estimator
    @unpack Nx = estimator.params.mesh
    @unpack sL, sR, indices, utilde = estimator.cache
    indices .= mod1.(j-2*sL-sR+2 : j+sL+2*sR-1, Nx)
    @views ushort = selectdim(uinit, 1, indices)
    if ndims(uinit)==1 # Scalar equations
        for i in 1:sL+sR-1
            utilde[i] = ushort[sL+sR]
        end
        for i in sL+sR:2*(sL+sR)-1
            utilde[i] = ushort[i]
        end
        for i in 2*(sL+sR):3*(sL+sR)-2
            utilde[i] = ushort[2*(sL+sR)-1]
        end
    # else    # Systems
    #     for r in 1:equation.p
    #         for i in 1:sL+sR-1
    #             utilde[i,r] = K[r]
    #         end
    #         for i in sL+sR:2*(sL+sR)-1
    #             utilde[i,r] = ushort[i,r]
    #         end
    #         for i in 2*(sL+sR):3*(sL+sR)-1
    #             utilde[i,r] = K[r]
    #         end
    #     end
    end
end