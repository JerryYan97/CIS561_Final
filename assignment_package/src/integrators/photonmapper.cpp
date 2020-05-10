#include "photonmapper.h"
#include "directlightingintegrator.h"
// #include <QMutex>
// QMutex mutex;
PhotonMapper::PhotonMapper(int numPhotons, std::vector<Photon> *photons, std::vector<Photon>* causticPhotons, Scene *s, std::shared_ptr<Sampler> sampler, int recursionLimit)
    : Integrator(Bounds2i(Point2i(0,0), Point2i(0,0)), s, sampler, recursionLimit), preprocessing(true), numPhotons(numPhotons), photons(photons), causticPhotons(causticPhotons)
{}

PhotonMapper::PhotonMapper(Bounds2i bounds, Scene *s, std::shared_ptr<Sampler> sampler, int recursionLimit, KDTree *kdtree, KDTree *causticKDtree, float searchR)
    : Integrator(bounds, s, sampler, recursionLimit), preprocessing(false), numPhotons(0), photons(nullptr), kdtree(kdtree), causticKDTree(causticKDtree), searchR(searchR)
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
        bool prevSpecBounce = false;
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

                            // change photon position
                            thisPhoton.pos = isect.point;

                            // flags variable shows whether this hit is a specular material.
                            // We still needs to know wheter previous hit is a specular material.
                            // store photon if non-specular and not first hit ?
                            // if ((flags & BSDF_SPECULAR) != 1 && j > 0)
                            // store photon in indirect photon map
                            if (!(flags & BSDF_SPECULAR))
                            {
                                thisPhoton.color *= (throughput / (float)scene->lights.size());
                                thisPhoton.wi = woW;
                                thisPhoton.hitPrim = isect.objectHit;
                                thisPhoton.hitNormal = isect.normalGeometric;

                                // If previous hit is a specular material, then this photon would be stored
                                // in the caustic photon tree.
                                if(prevSpecBounce)
                                {
                                    causticPhotons->push_back(thisPhoton);
                                }
                                else
                                {
                                    photons->push_back(thisPhoton);
                                }
                                // Owing to the fact that this hit is a diffuse material, we need to reset the spec bounce flag.
                                prevSpecBounce = false;
                            }
                            else
                            {
                                // We hit a specular material this time.
                                // We need to set the spec bounce flag.
                                // The throughput should only multiply albedo.
                                // Now, we don't care about it.
                                // f = isect.bsdf->Sample_f(wiW, &woW, sampler->Get2D(), &pdf, BSDF_ALL, &flags);
                                // throughput *= (f * AbsDot(woW, isect.bsdf->normal) / pdf);
                                // throughput *= (f / pdf);
                                if(MaxComponent(throughput) > 1.f || pdf == 0.f)
                                {
                                    int alart = 10;
                                }
                                prevSpecBounce = true;
                                ray = isect.SpawnRay(wiW);
                                continue;
                                // --j;
                            }
                            // change photon outgoing direction
                            // change photon radience
                            // change photon throughput
                            thisPhoton.wi = wiW;
                            f = isect.bsdf->f(wiW, woW);
                            pdf = isect.bsdf->Pdf(wiW, woW);
                            // float tempDot = AbsDot(woW, isect.bsdf->normal);
                            throughput *= (f * AbsDot(woW, isect.bsdf->normal) / pdf);
                            if(MaxComponent(throughput) > 1.f || pdf == 0.f)
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
    Color3f beta(1.f); // Throughput for specular material;
    if (scene.Intersect(ray, &isect))
    {
        if (isect.objectHit->material == nullptr)
        {
            result += isect.Le(woW);
        }
        else
        {
            isect.ProduceBSDF();
            float rangeArea = (1.f / 3.f) * Pi * searchR * searchR;
            // If hit a specular surface, then we are going to find another intersection point.
            // It can be light or other specular surface.
            // Therefore, we need to be careful.
            // Attention: we cannot reuse the 'isect' variable.
            // As long as current hit is a specular material:
            Ray mRay = ray;
            Vector3f wiW(0.f);
            float mPDF;
            BxDFType mFlag;
            Color3f mF;
            mF = isect.bsdf->Sample_f(woW, &wiW, sampler->Get2D(), &mPDF, BSDF_ALL, &mFlag);
            int specCounter = 0;
            std::vector<Vector3f> checkWiVec;
            checkWiVec.push_back(wiW);
            if(mFlag & BxDFType::BSDF_SPECULAR)
            {

                // mRay = isect.SpawnRay(wiW);
                mRay = Ray(isect.point + 0.1f * wiW, wiW);
                beta *= ((mF * AbsDot(woW, isect.bsdf->normal)) / mPDF);
                while(mFlag & BxDFType::BSDF_SPECULAR)
                {
                    specCounter++;
                    if(specCounter >= 2)
                    {
                        int alart = 10;
                    }

                    // beta *= ((mF * AbsDot(woW, mIsect.)));

                    Intersection mIsect;
                    woW = -wiW;


                    if(scene.Intersect(mRay, &mIsect))
                    {
                        if(mIsect.objectHit->material == nullptr)
                        {
                            return mIsect.Le(woW);
                        }
                    }
                    else
                    {
                        return Color3f(0.f);
                    }
                    // produce BSDF:
                    mIsect.ProduceBSDF();

                    mF = mIsect.bsdf->Sample_f(woW, &wiW, sampler->Get2D(), &mPDF, BSDF_ALL, &mFlag);

                    checkWiVec.push_back(wiW);
                    if(mFlag & BxDFType::BSDF_SPECULAR)
                    {
                        int alart = 10;
                        beta *= ((mF * AbsDot(woW, mIsect.bsdf->normal)) / mPDF);
                    }
                    if(!(mFlag & BxDFType::BSDF_SPECULAR))
                    {
                        isect = mIsect;
                        break;
                    }
                    mRay = mIsect.SpawnRay(wiW);
                }
            }
            /*
            mRay = isect.SpawnRay(wiW);
            while(mFlag & BxDFType::BSDF_SPECULAR)
            {
                Intersection mIsect;
                // Produce a new ray according to the sample_f:
                woW = -wiW;
                beta *= ((mF * AbsDot(wiW, isect.bsdf->normal)) / mPDF);
                // Find a new intersection.
                // If it is out of the scene, then return 0.f;
                // If it is a light source, then return light radiance;
                if(scene.Intersect(mRay, &mIsect))
                {
                    if(isect.objectHit->material == nullptr)
                    {
                        return isect.Le(woW);
                    }
                }
                else
                {
                    return Color3f(0.f);
                }
                // produce BSDF:
                mIsect.ProduceBSDF();
                mF = mIsect.bsdf->Sample_f(woW, &wiW, sampler->Get2D(), &mPDF, BSDF_ALL, &mFlag);
                if(!(mFlag & BxDFType::BSDF_SPECULAR))
                {
                    isect = mIsect;
                    break;
                }
            }
            */
            // beta = Color3f(1.f);

            // Direct light photons and indirect light photons rendering pass:
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
                    float disWeight = 1 - (glm::length(p->pos - isect.point) / searchR);
                    // float disWeight = 1.f;
                    Color3f thisPhoton = (pf * p->color * disWeight) * AbsDot(isect.normalGeometric, p->wi) / pdf;
                    photonResult += thisPhoton;
                    if(MaxComponent(thisPhoton) > 1.f || MinComponent(thisPhoton) < 0.f)
                    {
                        int alart = 10;
                    }
                }
            }
            if(countedPhoton != 0)
            {
                // photonResult = photonResult / (float)countedPhoton;
                // photonResult = photonResult;
                photonResult /= rangeArea;
                result = photonResult * beta;
                // result += (photonResult / rangeArea);
            }

            // Caustic photons rendering pass:
            std::vector<const Photon*> rangeCausticPhotons = causticKDTree->particlesInSphere(isect.point, searchR);
            Color3f causticPhotonResult(0.f);
            int countedCausticPhoton = 0;
            for(auto p : rangeCausticPhotons)
            {
                float pdf = isect.bsdf->Pdf(woW, p->wi);
                Color3f pf = isect.bsdf->f(woW, p->wi);
                if ((p->hitPrim == isect.objectHit) && Parallel(p->hitNormal, isect.normalGeometric))
                {
                    ++countedCausticPhoton;
                    if(pdf < 0.f)
                    {
                        continue;
                    }
                    if(IsBlack(pf) || pdf == 0.f)
                    {
                        continue;
                    }
                    float disWeight = 1 - (glm::length(p->pos - isect.point) / searchR);
                    Color3f thisPhoton = (pf * p->color * disWeight) * AbsDot(isect.normalGeometric, p->wi) / pdf;
                    causticPhotonResult += thisPhoton;
                    if(MaxComponent(thisPhoton) > 1.f || MinComponent(thisPhoton) < 0.f)
                    {
                        int alart = 10;
                    }
                }
            }
            if(countedCausticPhoton != 0)
            {
                causticPhotonResult /= rangeArea;
                result += (causticPhotonResult * beta);
                // result += Color3f(8.f);
            }
        }
    }
    return result;
}

