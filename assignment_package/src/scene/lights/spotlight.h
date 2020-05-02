#ifndef SPOTLIGHT_H
#define SPOTLIGHT_H
#include "light.h"

class SpotLight : public Light
{
public:
    SpotLight(const Transform &t, const Color3f& Le, float totalWidth, float falloffStart)
        :Light(t), emittedLight(Le),
          cosTotalWidth(std::cos(Radians(totalWidth))),
          cosFalloffStart(std::cos(Radians(falloffStart)))
    {}

    ~SpotLight(){}

    virtual Color3f L(const Intersection &isect, const Vector3f &w) const;

    Color3f Sample_Li(const Intersection &ref, const Point2f &xi,
                                         Vector3f *wi, Float *pdf) const;

    virtual float Pdf_Li(const Intersection &ref, const Vector3f &wi) const;

    float Falloff(const Vector3f &w) const;

    const Color3f emittedLight;
    const float cosTotalWidth, cosFalloffStart;
};

#endif // SPOTLIGHT_H
