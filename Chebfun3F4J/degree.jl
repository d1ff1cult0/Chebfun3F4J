# Heel vreemd dat hier U V V staat en niet U V W zoals ik zou denken, maar dit is ook hoe het in 
# de github van chebfun3F staat
function degree(cf3F)
    return [length(cf3F.U), length(cf3F.V), length(cf3F.V)]
end