#include "photonmapper.h"
#include "directlightingintegrator.h"

PhotonMapper::PhotonMapper(int numPhotons, std::vector<Photon> *photons, Scene *s, std::shared_ptr<Sampler> sampler, int recursionLimit)
    : Integrator(Bounds2i(Point2i(0,0), Point2i(0,0)), s, sampler, recursionLimit), preprocessing(true), numPhotons(numPhotons), photons(photons)
{}

PhotonMapper::PhotonMapper(Bounds2i bounds, Scene *s, std::shared_ptr<Sampler> sampler, int recursionLimit, KDTree *kdtree, float searchR)
    : Integrator(bounds, s, sampler, recursionLimit), preprocessing(false), numPhotons(0), photons(nullptr), kdtree(kdtree), searchR(searchR)
{}

void PhotonMapper::Render()
{
    /* Cast Photons */
    /* Use the idea in Naive Integrator */
    if(preprocessing)
    {
        // Determine how many photons to assign to each light source
        // given numPhotons and the intensity of each light.
        // Shoot a number of photons equal to numPhotons from
        // the lights, bouncing them through the scene and pushing
        // back the result of each bounce to the photons vector
        // stored in the PhotonMapper.
        int photonEachLight = numPhotons / scene->lights.size();
        for(auto l : scene->lights)
        {
            for(int i = 0; i < photonEachLight; ++i)
            {
                // 1. generate a set of photons for each light.
                int totalPhotons = 8 * photonEachLight;
                Photon thisPhoton = l->Sample_Photon(sampler, totalPhotons);
                Ray ray(thisPhoton.pos + 0.00001f * thisPhoton.wi, thisPhoton.wi);
                // thisPhoton.wi = - thisPhoton.wi;
                Color3f throughput(1.f);
                for(int j = 0; j < recursionLimit; ++j)
                {
                    // reverse the incoming photon direction
                    Vector3f woW = -glm::normalize(ray.direction);
                    Intersection isect;
                    if (scene->Intersect(ray, &isect))
                    {
                        if (isect.objectHit->material != nullptr)
                        {
                            isect.ProduceBSDF();
                            Vector3f wiW(0.f);
                            float pdf = 0.f;
                            BxDFType flags;
                            Color3f f = isect.bsdf->Sample_f(woW, &wiW, sampler->Get2D(), &pdf, BSDF_ALL, &flags);
                            wiW = glm::normalize(wiW);
                            if (IsBlack(f) || pdf == 0.f) break;
                            // Check whether hit a specular object.
                            // Then it is a caustic photon.
                            if(flags & BxDFType::BSDF_SPECULAR)
                            {

                            }


                            // change photon position
                            thisPhoton.pos = isect.point;

                            // store photon if non-specular and not first hit ?
                            // if ((flags & BSDF_SPECULAR) != 1 && j > 0)
                            // store photon in indirect photon map
                            if ((flags & BSDF_SPECULAR) != 1)
                            {
                                thisPhoton.color *= (throughput / (float)scene->lights.size());
                                thisPhoton.wi = woW;
                                thisPhoton.hitPrim = isect.objectHit;
                                thisPhoton.hitNormal = isect.normalGeometric;
                                photons->push_back(thisPhoton);
                            }
                            // change photon outgoing direction
                            // change photon radience
                            // change photon throughput
                            thisPhoton.wi = wiW;
                            f = isect.bsdf->f(wiW, woW);
                            pdf = isect.bsdf->Pdf(wiW, woW);
                            // float tempDot = AbsDot(woW, isect.bsdf->normal);
                            throughput *= (f * AbsDot(woW, isect.bsdf->normal) / pdf);
                            if(MaxComponent(throughput) > 1.f)
                            {
                                int alart = 10;
                            }
                            ray = isect.SpawnRay(wiW);

                            // possibly terminate the path with Russian roulette
                            /*if (j > 3)
                            {
                                float maxComp = MaxComponent(throughput);
                                float q = std::max(0.05f, 1.f - maxComp);
                                if (sampler->Get1D() < q) break;
                            }*/
                            // float maxComp = MaxComponent(throughput);
                            // float zeta = sampler->Get1D();
                            // if(maxComp < (1.f - zeta))
                            // {
                            //     break;
                            // }
                            // else
                            // {
                            //    throughput *= (1.f / maxComp);
                            // }
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
        }
    }
    else
    {
        Integrator::Render(); // Invokes Li for each ray from a pixel
    }
}

Color3f PhotonMapper::Li(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler, int depth) const
{
    Color3f result(0.f);
    Vector3f woW = -ray.direction;
    Intersection isect;
    if (scene.Intersect(ray, &isect))
    {
        if (isect.objectHit->material == nullptr)
        {
            result += isect.Le(woW);
        }
        else
        {
            isect.ProduceBSDF();
            float rangeArea = Pi * searchR * searchR;
            std::vector<const Photon*> rangePhotons = kdtree->particlesInSphere(isect.point, searchR);
            Color3f photonResult(0.f);
            int countedPhoton = 0;
            for (auto p : rangePhotons)
            {
                float pdf = isect.bsdf->Pdf(woW, p->wi);
                Color3f pf = isect.bsdf->f(woW, p->wi);
                if ((p->hitPrim == isect.objectHit) && Parallel(p->hitNormal, isect.normalGeometric))
                {
                    ++countedPhoton;
                    if(pdf < 0.f)
                    {
                        continue;
                    }
                    if(IsBlack(pf) || pdf == 0.f)
                    {
                        continue;
                    }
                    // float disWeight = 1 - (glm::length(p->pos - isect.point) / searchR);
                    float disWeight = 1.f;
                    Color3f thisPhoton = (pf * p->color * disWeight) * AbsDot(isect.normalGeometric, p->wi) / pdf;
                    photonResult += thisPhoton;
                    if(MaxComponent(thisPhoton) > 1.f)
                    {
                        int alart = 10;
                    }
                }
            }
            if(countedPhoton != 0)
            {
                photonResult = photonResult / (float)countedPhoton;
                // photonResult = photonResult;
                // photonResult /= rangeArea;
                result = photonResult;
                // result += (photonResult / rangeArea);
            }
        }
    }
    return result;
}

