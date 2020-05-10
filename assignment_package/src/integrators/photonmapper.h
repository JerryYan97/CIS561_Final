#pragma once
#include "integrator.h"
#include <scene/Photon.h>
#include <scene/kdtree.h>

class PhotonMapper : public Integrator
{
public:
    PhotonMapper(int numPhotons, std::vector<Photon>* photons, std::vector<Photon>* causticPhotons, Scene* s, std::shared_ptr<Sampler> sampler, int recursionLimit);
    PhotonMapper(Bounds2i bounds, Scene* s, std::shared_ptr<Sampler> sampler, int recursionLimit, KDTree *kdtree, KDTree *causticKDtree, float searchR);
    virtual void Render();
    virtual Color3f Li(const Ray& ray, const Scene& scene, std::shared_ptr<Sampler> sampler, int depth) const;
private:
    bool preprocessing;
    int numPhotons;
    std::vector<Photon>* photons;
    std::vector<Photon>* causticPhotons;

    float searchR;
    KDTree *kdtree;
    KDTree *causticKDTree;
};

