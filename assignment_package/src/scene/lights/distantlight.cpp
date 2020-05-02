#include "distantlight.h"

float DistantLight::Falloff(const Vector3f &w) const
{
    return glm::pow((worldRadius - glm::length(w - worldCenter)) / (worldRadius), 15.f);
}

Color3f DistantLight::L(const Intersection &isect, const Vector3f &w) const
{
    return Color3f(0.f);
}

Color3f DistantLight::Sample_Li(const Intersection &ref, const Point2f &xi, Vector3f *wi, Float *pdf) const
{
    *wi = -wLight;
    *pdf = 1;
    return this->Falloff(ref.point) * emittedLight;
}

float DistantLight::Pdf_Li(const Intersection &ref, const Vector3f &wi) const
{
    return 0;
}
