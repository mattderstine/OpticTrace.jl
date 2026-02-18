# Helper functions for testing

"""
define function to compute the normal vector using the gradient of the sag function
    This function takes a ray and a surface profile as input, and returns the normal vector at the intersection point of the ray with the surface.
The function first finds the intersection point of the ray with the surface using the deltaToSurf method, then finds the sag value at the intersection point, and finally finds the gradient of the sag value at the intersection point using an autodiff package. The normal vector is then obtained by normalizing the gradient vector.
"""
function normal_from_sag(ray::Ray, surfProfile)
    # find the intersection point of the ray with the surface using the deltaToSurf method
    delta = OpticTrace.deltaToSurf(ray, surfProfile)
    intersection_point = ray.base + delta * ray.dir
    # find the gradient of the sag value at the intersection point using an autodiff package
    grad_sag = ForwardDiff.gradient((x) -> sag(x[1], x[2], surfProfile), intersection_point[1:2])
    # the normal vector is then obtained by normalizing the gradient vector
    normal_vector = -normalize(Vec3(grad_sag[1], grad_sag[2], -1.0))
    return normal_vector
end

function normal_from_sag(intersection_point::Point3, surfProfile)
    # find the gradient of the sag value at the intersection point using an autodiff package
    grad_sag = ForwardDiff.gradient((x) -> sag(x[1], x[2], surfProfile), intersection_point[1:2])
    # the normal vector is then obtained by normalizing the gradient vector
    normal_vector = -normalize(Vec3(grad_sag[1], grad_sag[2], -1.0))
    return normal_vector
end

function normal_from_sag(x1::T, y1::T, surfProfile::S) where{T<:Real, S<:OpticTrace.AbstractSurfProfile{T}}
    f(x) = sag(x[1], x[2], surfProfile)

    x0 = [x1, y1]
    # find the gradient of the sag value at the intersection point using an autodiff package
    grad_sag = ForwardDiff.gradient(f, x0)
    # the normal vector is then obtained by normalizing the gradient vector
    normal_vector = -normalize(Vec3(grad_sag[1], grad_sag[2], -1.0))
    return normal_vector
end
function rprop(ray::Ray, delta::Float64)
    intersection_point = Point3(ray.base + delta * ray.dir)
    return intersection_point
end