#ifndef SPPMINTEGRATOR_H
#define SPPMINTEGRATOR_H
#include "integrator.h"

class SPPMIntegrator : public Integrator
{
public:
    SPPMIntegrator(Bounds2i bounds, Scene* s,
                   std::shared_ptr<Sampler> sampler, int nIterations,
                   int photonsPerIteration, int maxDepth, float initialSearchRadius, int writeFrequency)
        : Integrator(bounds, s, sampler, nIterations),
          initialSearchRadius(initialSearchRadius),
          maxDepth(maxDepth),
          photonsPerIteration(photonsPerIteration > 0
                              ? photonsPerIteration
                              : film->bounds.Area()),
          writeFrequency(writeFrequency),
          nIterations(nIterations)
    {}
    virtual void Render() override;
    virtual Color3f Li(const Ray& ray, const Scene& scene, std::shared_ptr<Sampler> sampler, int depth) const override
    {
        return Color3f(0.f, 0.f, 0.f);
    }
    ~SPPMIntegrator(){}
private:
    const Float initialSearchRadius;
    const int maxDepth;
    const int photonsPerIteration;
    const int writeFrequency;
    const int nIterations;
};

#endif // SPPMINTEGRATOR_H
