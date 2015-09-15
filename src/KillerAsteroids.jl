module KillerAsteroids

require(joinpath(Pkg.dir("KillerAsteroids"), "src", "WisePSF.jl"))
require(joinpath(Pkg.dir("KillerAsteroids"), "src", "ModelInit.jl"))
require(joinpath(Pkg.dir("KillerAsteroids"), "src", "Synthetic.jl"))
require(joinpath(Pkg.dir("KillerAsteroids"), "src", "Probability.jl"))

import WisePSF
import Synthetic
import Probability

using ModelInit

export Image, AsteroidParams, Prio

end # module
