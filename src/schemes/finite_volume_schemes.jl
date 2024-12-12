struct FVScheme <: Scheme
    timeScheme::TimeScheme
    spaceScheme::SpaceScheme
end

get_name(fvs::FVScheme) = get_name(fvs.timeScheme)*" "*get_name(fvs.spaceScheme)

get_sL(fvs::FVScheme) = get_sL(fvs.timeScheme, fvs.spaceScheme)
get_sR(fvs::FVScheme) = get_sR(fvs.timeScheme, fvs.spaceScheme)