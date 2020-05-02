#include "fulllightingintegrator.h"
#include "directlightingintegrator.h"
#include <scene/lights/pointlight.h>

Color3f FullLightingIntegrator::Li(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler, int depth) const
{
    // Instantiate an accumlated ray color that begins as black;
    Color3f L(0.f);
    // Instantiate an accumulated ray throughput color that begins as white;
    // The throughput will be used to determine when your ray path terminates via the Russian Roulette heuristic.
    Color3f beta(1.f);

    Ray rayPath(ray);

    bool specBounce = false;

    int bounceCounter = depth;

    while(bounceCounter > 0)
    {
        Intersection isect;
        if(!scene.Intersect(rayPath, &isect))
        {
            break;
        }

        // Initialize common variable for integrator.
        Vector3f wo = - rayPath.direction;

        // Compute emitted light if ray hit an area light source.
        // if(isect.objectHit->GetAreaLight() || (specBounce && (depth - bounceCounter == 1)))
        if(isect.objectHit->GetAreaLight())
        {
            if(bounceCounter == depth || specBounce)
            {
                Color3f lightSource = isect.Le(wo);
                L += beta * lightSource;
                break;
            }
            else
            {
                break;
            }
        }
        // Ask _objectHit_ to produce a BSDF
        isect.ProduceBSDF();

        // If previous hit is a specular point, while this hit is not light source.
        // Then, we reset the specBounce flag.
        specBounce = false;

        // Initialize normal for integrator:
        // Normal3f n = isect.normalGeometric;
        Normal3f n = isect.bsdf->normal;

        // Check whether hit a specular object
        Color3f mLiSpec(0.f);
        Vector3f wiSpec;
        float pdfSpec;
        BxDFType flagsSpec;
        Color3f fSpec = isect.bsdf->Sample_f(wo, &wiSpec, sampler->Get2D(), &pdfSpec, BSDF_ALL, &flagsSpec);
        if(flagsSpec & BxDFType::BSDF_SPECULAR)
        {
            beta *= (fSpec * AbsDot(wiSpec, n) / pdfSpec);
            specBounce = true;
            rayPath = isect.SpawnRay(wiSpec);
            --bounceCounter;
            continue;
        }

        Color3f LTerm = DirectLight(rayPath, scene, sampler);
        L += (beta * LTerm);

        // Computing the ray bounce and global illumination.
        Color3f mLiG(0.f);
        Vector3f wiG;
        float pdfG;
        BxDFType flagsG;
        Color3f fG = isect.bsdf->Sample_f(wo, &wiG, sampler->Get2D(), &pdfG, BSDF_ALL, &flagsG);
        if(IsBlack(fG) || pdfG == 0.f)
        {
            break;
        }
        beta *= fG * AbsDot(wiG, n) / pdfG;
        rayPath = isect.SpawnRay(wiG);

        // Correctly accounting for direct lighting.

        // Russian Roulette Ray Termination.
        // Compare the maximum RGB component of your throughput to a uniform random number and
        // stop your while loop if said component is smaller than the random number.
        float maxChannel = beta[0];
        for(int i = 1; i < 3; i++)
        {
            if(beta[i] > maxChannel)
            {
                maxChannel = beta[i];
            }
        }

        float zeta = sampler->Get1D();
        if(maxChannel < (1.f - zeta))
        {
            break;
        }
        else
        {
            beta *= (1.f / maxChannel);
        }

        --bounceCounter;
    }



    return L;
}
