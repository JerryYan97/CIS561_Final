#include "directlightingintegrator.h"
#include <scene/lights/pointlight.h>
#include <scene/lights/spotlight.h>
#include <scene/lights/distantlight.h>
#include <iostream>

float BalanceHeuristic(int nf, Float fPdf, int ng, Float gPdf)
{
    //TODO
    return (float(nf) * fPdf) / (float(nf) * fPdf + float(ng) * gPdf);
}

float PowerHeuristic(int nf, Float fPdf, int ng, Float gPdf)
{
    float f = nf * fPdf;
    float g = ng * gPdf;
    //TODO
    return (f * f) / (f * f + g * g);
}

Color3f DirectLightingIntegrator::Li(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler, int depth) const
{
    //TODO
    Color3f L(0.f);


    // Find closest ray intersection or return background radiance.
    Intersection isect;

    if(!scene.Intersect(ray, &isect))
    {
        return L;
    }

    // Compute emitted and reflected light at ray intersection point

    // Initialize common variable for Whitted integrator
    Normal3f n = isect.normalGeometric;
    Vector3f wo = - ray.direction;

    // Compute scattering functions for surface interaction

    // Compute emitted light if ray hit an area light source
    Color3f lightSource = isect.Le(wo);
    if(!IsBlack(lightSource))
    {
        L += lightSource;
        return L;
    }

    // Ask _objectHit_ to produce a BSDF
    // based on other data stored in this
    // Intersection, e.g. the surface normal
    // and UV coordinates
    isect.ProduceBSDF();

    // Light PDF sampling:
    // Randomly select a light source from scene.lights and call its Sample_Li function:
    int nLights = int(scene.lights.size());
    if (nLights == 0) return Color3f(0.f);
    int lightNum = std::min((int)(sampler->Get1D() * nLights), nLights - 1);
    const std::shared_ptr<Light> &light = scene.lights[lightNum];

    Vector3f wi(0.f);
    float pdf = 0.f;
    Color3f Li = light->Sample_Li(isect, sampler->Get2D(), &wi, &pdf);

    // Determine the type of light source:
    Light* lightPtr = light.get();
    PointLight* pLightPtr = dynamic_cast<PointLight*>(lightPtr);
    SpotLight* spLightPtr = dynamic_cast<SpotLight*>(lightPtr);
    DistantLight* disLightPtr = dynamic_cast<DistantLight*>(lightPtr);
    if(pLightPtr || spLightPtr || disLightPtr)
    {
        Ray tempRay(Vector3f(0, 0, 0), Vector3f(0, 0, 0));

        if(pLightPtr || spLightPtr)
        {
            // Shadow Test from point light to the reference point:
            tempRay = Ray(lightPtr->transform.position(), -wi);
        }
        else if(disLightPtr)
        {
            // Shadow Test from a point 2 radius away to the reference point:
            tempRay = Ray(isect.point + 2.f * wi * disLightPtr->worldRadius, -wi);
        }

        Intersection tempInsect;
        if(!scene.Intersect(tempRay, &tempInsect))
        {
            Li = Color3f(0.f);
        }
        else
        {
            Color3f f = isect.bsdf->f(wo, wi);
            float lMaterialPDF = isect.bsdf->Pdf(wo, wi);
            if(pdf <= 0.f || lMaterialPDF <= 0.f)
            {
                Li = Color3f(0.f);
            }
            else
            {
                if(glm::length(tempInsect.point - isect.point) < 0.01)
                {
                    // Li = glm::clamp(Li, 0.f, 1.f);
                    L += (f * AbsDot(wi, n) * Li / (pdf / nLights));
                }
            }
        }
        return L;
    }
    else
    {
        // Shadow Test and evaluate the LTE for Light PDF sampling:
        Ray tempRay = Ray(isect.point + isect.normalGeometric * 1e-4f, wi);
        Intersection tempInsect;
        if(!scene.Intersect(tempRay, &tempInsect))
        {
            Li = Color3f(0.f);
        }
        else
        {
            if(light.get() != tempInsect.objectHit->GetAreaLight())
            {
                Li = Color3f(0.f);
            }
            else
            {
                Color3f f = isect.bsdf->f(wo, wi);
                float lMaterialPDF = isect.bsdf->Pdf(wo, wi);
                if(pdf <= 0.f || lMaterialPDF <= 0.f)
                {
                    Li = Color3f(0.f);
                }
                else
                {
                    // W0 : Wj = 1 : nf, ng.
                    float lightPDFWeight = PowerHeuristic(1, pdf, 1, lMaterialPDF);
                    // float lightPDFWeight = 1.f;
                    // std::cout << "lightPDFWeight:" << lightPDFWeight << std::endl;
                    L += f * Li * AbsDot(wi, n) * lightPDFWeight / (pdf / nLights);
                }
            }
        }
        // ************************************************* //

        // BRDF PDF sampling:
        Vector3f mWiB(0.f);
        Color3f mLi(0.f);
        float mPDFMaterial = 0.f;
        BxDFType flags;
        Color3f mF = isect.bsdf->Sample_f(wo, &mWiB, sampler->Get2D(), &mPDFMaterial, BSDF_ALL, &flags);
        if(!IsBlack(mF) && mPDFMaterial != 0.f)
        {
            Ray tempRay = isect.SpawnRay(mWiB);
            Intersection insectBRDF;
            if(scene.Intersect(tempRay, &insectBRDF))
            {
                if(insectBRDF.objectHit->GetAreaLight() == light.get())
                {
                    const AreaLight* arealightPtr = insectBRDF.objectHit->GetAreaLight();
                    float pdfLWB = arealightPtr->Pdf_Li(isect, mWiB);
                    float wWiB = PowerHeuristic(1, mPDFMaterial, 1, pdfLWB);
                    mLi = mF * Li * AbsDot(mWiB, n) * wWiB / (pdfLWB / nLights);
                    L += mLi;
                }
            }
        }

        if(MaxComponent(L) > 1.f)
        {
            int alart = 10;
        }

        return L;
    }
}


Color3f DirectLight(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler)
{
    //TODO
    Color3f L(0.f);


    // Find closest ray intersection or return background radiance.
    Intersection isect;

    if(!scene.Intersect(ray, &isect))
    {
        return L;
    }

    // Compute emitted and reflected light at ray intersection point

    // Initialize common variable for Whitted integrator
    Normal3f n = isect.normalGeometric;
    Vector3f wo = - ray.direction;

    // Compute scattering functions for surface interaction

    // Compute emitted light if ray hit an area light source
    Color3f lightSource = isect.Le(wo);
    if(!IsBlack(lightSource))
    {
        L += lightSource;
        return L;
    }

    // Ask _objectHit_ to produce a BSDF
    // based on other data stored in this
    // Intersection, e.g. the surface normal
    // and UV coordinates
    isect.ProduceBSDF();

    // Light PDF sampling:
    // Randomly select a light source from scene.lights and call its Sample_Li function:
    int nLights = int(scene.lights.size());
    if (nLights == 0) return Color3f(0.f);
    int lightNum = std::min((int)(sampler->Get1D() * nLights), nLights - 1);
    const std::shared_ptr<Light> &light = scene.lights[lightNum];

    Vector3f wi(0.f);
    float pdf = 0.f;
    Color3f Li = light->Sample_Li(isect, sampler->Get2D(), &wi, &pdf);

    // Determine the type of light source:
    Light* lightPtr = light.get();
    PointLight* pLightPtr = dynamic_cast<PointLight*>(lightPtr);
    SpotLight* spLightPtr = dynamic_cast<SpotLight*>(lightPtr);
    DistantLight* disLightPtr = dynamic_cast<DistantLight*>(lightPtr);
    if(pLightPtr || spLightPtr || disLightPtr)
    {
        // Shadow Test from point light to the reference point:
        Ray tempRay(Vector3f(0, 0, 0), Vector3f(0, 0, 0));

        if(pLightPtr || spLightPtr)
        {
            // Shadow Test from point light to the reference point:
            tempRay = Ray(lightPtr->transform.position(), -wi);
        }
        else if(disLightPtr)
        {
            // Shadow Test from a point 2 radius away to the reference point:
            tempRay = Ray(isect.point + 2.f * wi * disLightPtr->worldRadius, -wi);
        }

        Intersection tempInsect;
        if(!scene.Intersect(tempRay, &tempInsect))
        {
            Li = Color3f(0.f);
        }
        else
        {
            Color3f f = isect.bsdf->f(wo, wi);
            float lMaterialPDF = isect.bsdf->Pdf(wo, wi);
            if(pdf <= 0.f || lMaterialPDF <= 0.f)
            {
                Li = Color3f(0.f);
            }
            else
            {
                if(glm::length(tempInsect.point - isect.point) < 0.01)
                {
                    // Li = glm::clamp(Li, 0.f, 1.f);
                    L += (f * AbsDot(wi, n) * Li / (pdf / nLights));
                }
            }
        }
        return L;
        /*Ray tempRay = Ray(lightPtr->transform.position(), -wi);
        Intersection tempInsect;
        if(!scene.Intersect(tempRay, &tempInsect))
        {
            Li = Color3f(0.f);
        }
        else
        {
            if(glm::length(tempInsect.point - isect.point) < 0.01)
            {
                Li = glm::clamp(Li, 0.f, 1.f);
                L += (AbsDot(wi, n) * Li / float(nLights));
            }
        }
        return L;*/
    }
    else
    {
        // Shadow Test and evaluate the LTE for Light PDF sampling:
        Ray tempRay = Ray(isect.point + isect.normalGeometric * 1e-4f, wi);
        Intersection tempInsect;
        if(!scene.Intersect(tempRay, &tempInsect))
        {
            Li = Color3f(0.f);
        }
        else
        {
            if(light.get() != tempInsect.objectHit->GetAreaLight())
            {
                Li = Color3f(0.f);
            }
            else
            {
                Color3f f = isect.bsdf->f(wo, wi);
                float lMaterialPDF = isect.bsdf->Pdf(wo, wi);
                if(pdf <= 0.f || lMaterialPDF <= 0.f)
                {
                    Li = Color3f(0.f);
                }
                else
                {
                    // W0 : Wj = 1 : nf, ng.
                    float lightPDFWeight = PowerHeuristic(1, pdf, 1, lMaterialPDF);
                    // float lightPDFWeight = 1.f;
                    // std::cout << "lightPDFWeight:" << lightPDFWeight << std::endl;
                    L += f * Li * AbsDot(wi, n) * lightPDFWeight / (pdf / nLights);
                }
            }
        }
    }

    // ************************************************* //

    // BRDF PDF sampling:
    Vector3f mWiB(0.f);
    Color3f mLi(0.f);
    float mPDFMaterial = 0.f;
    BxDFType flags;
    Color3f mF = isect.bsdf->Sample_f(wo, &mWiB, sampler->Get2D(), &mPDFMaterial, BSDF_ALL, &flags);
    if(!IsBlack(mF) && mPDFMaterial != 0.f)
    {
        Ray tempRay = isect.SpawnRay(mWiB);
        Intersection insectBRDF;
        if(scene.Intersect(tempRay, &insectBRDF))
        {
            if(insectBRDF.objectHit->GetAreaLight() == light.get())
            {
                const AreaLight* arealightPtr = insectBRDF.objectHit->GetAreaLight();
                float pdfLWB = arealightPtr->Pdf_Li(isect, mWiB);
                float wWiB = PowerHeuristic(1, mPDFMaterial, 1, pdfLWB);
                mLi = mF * Li * AbsDot(mWiB, n) * wWiB / (pdfLWB / nLights);
                L += mLi;
            }
        }
    }

    return L;
}
