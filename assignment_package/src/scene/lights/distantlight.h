#ifndef DISTANTLIGHT_H
#define DISTANTLIGHT_H
#include "light.h"

class DistantLight  : public Light
{
public:
    DistantLight(const Transform &t, const Color3f& Le, const Vector3f &wLight, const float iRadius)
        : Light(t), emittedLight(Le), wLight(glm::normalize(t.T3() * wLight)), worldCenter(t.position()), worldRadius(iRadius)
    {}

    virtual Color3f L(const Intersection &isect, const Vector3f &w) const;

    Color3f Sample_Li(const Intersection &ref, const Point2f &xi,
                                         Vector3f *wi, Float *pdf) const;

    virtual float Pdf_Li(const Intersection &ref, const Vector3f &wi) const;

    float Falloff(const Vector3f &w) const;

    ~DistantLight(){}

    const Color3f emittedLight;
    const Vector3f wLight;
    Point3f worldCenter;
    float worldRadius;
};

#endif // DISTANTLIGHT_H
