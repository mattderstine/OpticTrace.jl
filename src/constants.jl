
const ORIGIN = SVector(0., 0., 0.)
const ZAXIS = SVector(0., 0., 1.)
const YAXIS = SVector(0.,1.,0.)
const XAXIS = SVector(1.,0.,0.)

const ∞ = Inf

global refIndexDefault::Float64 = 1.0

rInDef()=refIndexDefault
function setRInDef(in)
    global refIndexDefault = in
end

export ∞, ORIGIN, ZAXIS, YAXIS, XAXIS, refIndexDefault, rInDef, setRInDef