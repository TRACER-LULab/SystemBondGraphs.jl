using NetworkDynamics
using LightGraphs

function R(v, edges, p, t)
    e = v[1]
    f = v[2]
    R = p[1]
    e = R * f    
end

function C(dv, v, edges, p, t)
    e = v[1]
    f = v[2]
    q = v[3]
    de = dv[1]
    df = dv[2]
    dq = dv[3]
    dq = f
    q = C * e

end