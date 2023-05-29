module FastSSC

using LoopVectorization

function SSC_integral(volume_element, WA, WB, WC, WD, responseAB, responseCD, σ², nz, nl, dz)
    ntomo = 10
    result = zeros(length(ℓ_array), length(ℓ_array), ntomo, ntomo, ntomo, ntomo)

    @tturbo for idxli in 1:nl
        for idxlj in 1:nl
            for i in 1:ntomo
                for j in 1:ntomo
                    for k in 1:ntomo
                        for l in 1:ntomo
                            for idxzi in 1:nz
                                for idxzj in 1:nz
                                    result[idxli, idxlj, i, j, k, l] += volume_element[idxzi] * volume_element[idxzj] * WA[i, idxzi] * WB[j, idxzi] * WC[k, idxzj] * WD[l, idxzj] *
                                                                        responseAB[idxli, idxzi] * responseCD[idxlj, idxzj] * σ²[idxzi, idxzj]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    return (dz^2) .* result
end
end # module FastSSC
