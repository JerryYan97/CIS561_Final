#include "spotlight.h"

float SpotLight::Falloff(const Vector3f &w) const
{
    Vector3f wl = glm::normalize(transform.invT3() * w);
    float cosTheta = wl.z;
    if (cosTheta < cosTotalWidth) return 0;
    if (cosTheta >= cosFalloffStart) return 1;
    float delta = (cosTheta - cosTotalWidth) / (cosFalloffStart - cosTotalWidth);
    return (delta * delta) * (delta * delta);
}

Color3f SpotLight::L(const Intersection &isect, const Vector3f &w) const
{
    return Color3f(0.f);
}

Color3f SpotLight::Sample_Li(const Intersection &ref, const Point2f &xi,
                                     Vector3f *wi, Float *pdf) const
{
    *wi = glm::normalize(transform.position() - ref.point);
    *pdf = 1.f;
    return emittedLight * Falloff(-*wi) / DistanceSquared(transform.position(), ref.point);
}

float SpotLight::Pdf_Li(const Intersection &ref, const Vector3f &wi) const
{
    return 0;
}
