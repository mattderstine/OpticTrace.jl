
export ∞, ORIGIN, ZAXIS, YAXIS, XAXIS, refIndexDefault, rInDef, setRInDef

const ORIGIN = Point3(0., 0., 0.)
const ZAXIS = Vec3(0., 0., 1.)
const YAXIS = Vec3(0.,1.,0.)
const XAXIS = Vec3(1.,0.,0.)

const ∞ = Inf

global refIndexDefault::Float64 = 1.0

rInDef()=refIndexDefault
rInDef(wave) = refIndexDefault #used as a glass function in the glass dictionary
function setRInDef(in)
    global refIndexDefault = in
end



const identityPol = SMatrix{3, 3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
