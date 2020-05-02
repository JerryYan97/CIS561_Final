#include "sppmintegrator.h"
#include "directlightingintegrator.h"
#include <scene/lights/pointlight.h>
#include <iostream>

struct SPPMPixel{

    float radius = 0;
    bool hitLight = false;
    Point2i pixel;
    Color3f Ld = Color3f(0.f);
    int M = 0;
    Color3f Phi = Color3f(0.f);
    struct VisiblePoint{
        Intersection vpIsect;
        // Point3f p;
        Vector3f wo;
        // BSDF bsdf = BSDF(Intersection(), 1);
        Color3f beta = Color3f(0.f);
        bool exist = false;
    } vp;
    float N = 0;
    Color3f tau = Color3f(0.f);
};

struct SPPMPixelListNode{
    SPPMPixel *pixel;
    SPPMPixelListNode *next;
};




static bool ToGrid(const Point3f &p, const Bounds3f &bounds, const int gridRes[3], Point3i *pi)
{
    bool inBounds = true;
    Vector3f pg = bounds.Offset(p);
    for (int i = 0; i < 3; ++i) {
        (*pi)[i] = (int)(gridRes[i] * pg[i]);
        inBounds &= ((*pi)[i] >= 0 && (*pi)[i] < gridRes[i]);
        (*pi)[i] = glm::clamp((*pi)[i], 0, gridRes[i] - 1);
    }
    return inBounds;
}

inline unsigned int hash(const Point3i &p, int hashSize)
{
    return (unsigned int)((p.x * 73856093) ^ (p.y * 19349663) ^ (p.z * 83492791)) % hashSize;
}

void SPPMIntegrator::Render()
{
    // Initialize pixelBounds and pixels array for SPPM.
    std::vector<Point2i> tilePixels = bounds.GetPoints();
    int nPixels = bounds.Area();
    std::unique_ptr<SPPMPixel[]> pixels(new SPPMPixel[nPixels]);
    for(int i = 0; i < nPixels; ++i)
    {
        pixels[i].radius = initialSearchRadius;
    }


    // Compute lightDistr for sampling lights proportional to power.

    // Perform nIterations of SPPM integration.
    for(int iter = 0; iter < nIterations; ++iter)
    {
        /* Generate SPPM visible points */
        // Pixel idx:
        int pIdx = 0;
        for(Point2i pixel : tilePixels)
        {
            // Get SPPMPixel for current pixel:
            SPPMPixel &currSPPMPixel = pixels[pIdx];
            currSPPMPixel.pixel = pixel;

            // Generate camera ray for pixel for SPPM:
            Point2f sample = sampler->Get2D() + Point2f(pixel);
            Ray ray = camera->Raycast(sample);
            Color3f beta(1.f);

            // Follow camera ray path until a visible point is created:
            bool specularBounce = false;
            for(int depth = 0; depth < 10; ++depth)
            {
                Intersection isect;
                if(!scene->Intersect(ray, &isect))
                {
                    // Accumulate light contributions for ray with no intersection:
                    break;
                }

                Vector3f wo = -ray.direction;

                // Process SPPM camera ray intersection:
                // Compute BSDF at SPPM camera ray intersection:
                // if(!isect.ProduceBSDF())
                if(isect.objectHit->GetAreaLight())
                {
                    // Hit the area light:
                    currSPPMPixel.hitLight = true;
                    currSPPMPixel.Ld = beta * isect.Le(wo);
                    currSPPMPixel.vp = {isect, wo, beta, true};
                    break;
                }
                isect.ProduceBSDF();
                const BSDF &bsdf = *isect.bsdf;

                // Possibly create visible point and end camera path:
                // ? This is not feasible for this project.
                // ? It would leave visible points on glass material.
                // bool isDiffuse = bsdf.BxDFsMatchingFlags(BxDFType(BSDF_DIFFUSE | BSDF_REFLECTION | BSDF_TRANSMISSION)) > 0;
                // bool isGlossy = bsdf.BxDFsMatchingFlags(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)) > 0;
                // bool isDiffuse = bsdf.BxDFsMatchingFlags(BxDFType(BSDF_DIFFUSE | BSDF_GLOSSY)) > 0;
                // bool isGlossy = bsdf.BxDFsMatchingFlags(BxDFType(BSDF_GLOSSY));
                // if(isDiffuse)
                // {
                //
                // }

                // Spwan ray from SPPM camera path vertex
                // ? Spwan ray from specular material? or Transmissive material ?
                // ? There are only two cases for materials: diffuse or specular(include trans) ?
                float pdf;
                Vector3f wi;
                BxDFType type;
                Color3f f = bsdf.Sample_f(wo, &wi, sampler->Get2D(), &pdf, BSDF_ALL, &type);

                if(type & BxDFType::BSDF_SPECULAR)
                {
                    beta *= (f * AbsDot(wi, bsdf.normal) / pdf);
                    if(beta.y < 0.25)
                    {
                        float continueProb = std::min(1.f, beta.y);
                        if(sampler->Get1D() > continueProb)
                        {
                            break;
                        }
                        beta /= continueProb;
                    }
                    ray = isect.SpawnRay(wi);
                }
                else
                {
                    // Accumulate direct illumination at SPPM camera ray intersection:
                    currSPPMPixel.Ld += (beta * DirectLight(ray, *scene, sampler));
                    currSPPMPixel.vp = {isect, wo, beta, true};
                    // currSPPMPixel.vp = {isect.point, wo, &bsdf, beta, true};
                    break;
                }
            }
            ++pIdx;
        }

        // Create grid of all SPPM visible points:
        int hashSize = nPixels;
        std::vector<SPPMPixelListNode*> grid(hashSize);
        Bounds3f gridBounds;
        float maxRadius = 0.f;
        for(int i = 0; i < nPixels; ++i)
        {
            const SPPMPixel &currSPPMPixel = pixels[i];
            if(IsBlack(currSPPMPixel.vp.beta))
            {
                continue;
            }
            Bounds3f vpBound = Expand(Bounds3f(currSPPMPixel.vp.vpIsect.point), currSPPMPixel.radius);
            gridBounds = Union(gridBounds, vpBound);
            maxRadius = std::max(maxRadius, currSPPMPixel.radius);
        }
        // Compute resolution of SPPM grid in each dimension:
        Vector3f diag = gridBounds.Diagonal();
        float maxDiag = MaxComponent(diag);
        int baseGridRes = (int)(maxDiag / maxRadius);
        int gridRes[3];
        for(int i = 0; i < 3; ++i)
        {
            gridRes[i] = std::max((int)(baseGridRes * diag[i] / maxDiag), 1);
        }

        int validPointCounter = 0;
        // Add visible points to SPPM grid:
        for(int pixelIdx = 0; pixelIdx < nPixels; ++pixelIdx)
        {
            SPPMPixel& currSPPMPixel = pixels[pixelIdx];
            if(!IsBlack(currSPPMPixel.vp.beta))
            {
                // Add pixel's visible point to applicable grid cells:
                float radius = currSPPMPixel.radius;
                Point3i pMin, pMax;
                ToGrid(currSPPMPixel.vp.vpIsect.point - Vector3f(radius, radius, radius), gridBounds, gridRes, &pMin);
                ToGrid(currSPPMPixel.vp.vpIsect.point + Vector3f(radius, radius, radius), gridBounds, gridRes, &pMax);
                for(int z = pMin.z; z <= pMax.z; ++z)
                {
                    for(int y = pMin.y; y <= pMax.y; ++y)
                    {
                        for(int x = pMin.x; x <= pMax.x; ++x)
                        {
                            // Add visible point to grid cell (x, y, z)
                            int h = hash(Point3i(x, y, z), hashSize);
                            SPPMPixelListNode *node = new SPPMPixelListNode();
                            node->pixel = &currSPPMPixel;
                            // Atomically add node to the start of grid[h]'s linked list:
                            node->next = grid[h];
                            grid[h] = node;
                            ++validPointCounter;
                        }
                    }
                }
            }
        }

        /* Insanity Check */
        /* Print out every pixels in grid:
         * for(int idx = 0; idx < nPixels; ++idx)
        {
            if(grid[idx])
            {
                for(SPPMPixelListNode* ptr = grid[idx]; ptr != nullptr; ptr = ptr->next)
                {
                    film->SetPixelColor(ptr->pixel->pixel, glm::clamp(ptr->pixel->Ld, 0.f, 1.f));
                }
            }
        }
        */

        // Trace photons and accumulate contributions
        for(int photonIdx = 0; photonIdx < photonsPerIteration; ++photonIdx)
        {
            // uint64_t haltonIndex = (uint64_t)iter * (uint64_t)photonsPerIteration + photonIdx;
            // int haltonDim = 0;
            // Choose light to shoot photon from:
            // Follow Direct light integrator. Randomly select a light source from scene.
            int nLights = int(scene->lights.size());
            if (nLights == 0) return;
            int lightNum = std::min((int)(sampler->Get1D() * nLights), nLights - 1);
            const std::shared_ptr<Light> &light = scene->lights[lightNum];

            // Compute sample values for photon ray leaving light source:
            float lightPdf = 1.f / (4 * Pi); // Only for point light.

            // Generate photonRay from light source and initialize beta:
            Ray photonRay(Point3f(0.f), Vector3f(0.f));
            Normal3f nLight;
            float pdfPos, pdfDir;
            Color3f Le = light->Sample_Le(sampler->Get2D(), sampler->Get2D(), &photonRay, &nLight, &pdfPos, &pdfDir);
            if (pdfPos == 0 || pdfDir == 0 || IsBlack(Le)) return;
            Color3f beta = (AbsDot(nLight, photonRay.direction) * Le) / (lightPdf * pdfPos * pdfDir);
            if(IsBlack(beta))
            {
                return;
            }

            // Follow photon path through scene and record intersections:
            Intersection isect;
            for(int depth = 0; depth < 10; ++depth)
            {
                if(!scene->Intersect(photonRay, &isect))
                {
                    break;
                }
                // if(depth > 0)
                if(depth > 0)
                {
                    // Add photon contribution to nearby visible points:
                    // The photon makes no contribution at the first intersection.
                    Point3i photonGridIndex;
                    if(ToGrid(isect.point, gridBounds, gridRes, &photonGridIndex))
                    {
                        int h = hash(photonGridIndex, hashSize);
                        // Add photon contribution to visible points in grid[h].
                        for(SPPMPixelListNode *node = grid[h]; node != nullptr; node = node->next)
                        {
                            SPPMPixel &pixel = *node->pixel;
                            float radius = pixel.radius;
                            if(DistanceSquared(pixel.vp.vpIsect.point, isect.point) > radius * radius)
                            {
                                continue;
                            }
                            // Update pixel Phi and M for nearby photon.
                            Vector3f wi = -photonRay.direction;
                            // Color3f Phi = beta * pixel.vp.bsdf.f(pixel.vp.wo, wi);
                            if(pixel.vp.vpIsect.bsdf->numBxDFs == 0)
                            {
                                // Temporarily don't consider the area light without BSDF.
                                pixel.vp.vpIsect.ProduceBSDF();
                            }

                            Color3f Phi = beta * pixel.vp.vpIsect.bsdf->f(pixel.vp.wo, wi);
                            pixel.Phi += Phi;
                            ++pixel.M;
                        }
                    }
                }
                // Sample new photon ray direction:
                // Compute BSDF at photon intersection point:
                isect.ProduceBSDF();
                if(!isect.bsdf)
                {
                    --depth;
                    photonRay = isect.SpawnRay(photonRay.direction);
                    continue;
                }
                const BSDF &photonBSDF = *isect.bsdf;

                // Sample BSDF fr and direction wi for reflected photon:
                Vector3f wi, wo = -photonRay.direction;
                float pdf;
                BxDFType flags;
                // Generate bsdfSample for outgoing photon sample: ?
                Color3f fr = photonBSDF.Sample_f(wo, &wi, sampler->Get2D(), &pdf, BSDF_ALL, &flags);
                if(IsBlack(fr) || pdf == 0.f)
                {
                    break;
                }
                // Color3f bnew = beta * fr * AbsDot(wi, photonBSDF.normal) / pdf;
                beta = beta * fr * AbsDot(wi, photonBSDF.normal) / pdf;

                // Possibly terminate photon path with Russian Roulette:
                float maxChannel = beta[0];
                for(int i = 1; i < 3; ++i)
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
                    // float q = std::max(0.f, 1 - bnew.y / beta.y);
                    // beta = bnew / (1 - q);
                    beta *= (1.f / maxChannel);
                    photonRay = isect.SpawnRay(wi);
                }
            }
        }

        // Update pixel values from this pass's photons:
        for(int i = 0; i < nPixels; ++i)
        {
            SPPMPixel &p = pixels[i];
            if(p.M > 0)
            {
                // Update pixel photon count, search radius, and t from photons
                float gamma = 2.f / 3.f;
                float Nnew = p.N + gamma * p.M;
                float Rnew = p.radius * std::sqrt(Nnew / (p.N + p.M));
                Color3f Phi;
                Phi = p.Phi;
                p.tau = (p.tau + p.vp.beta * Phi) * (Rnew * Rnew) / (p.radius * p.radius);
                p.N = Nnew;
                p.radius = Rnew;

                p.M = 0;
                p.Phi = Color3f(0.f);
            }
            // Reset VisiblePoint in pixel
            p.vp.beta = Color3f(0.f);
            p.vp.vpIsect = Intersection();
            // p.vp.bsdf = BSDF(Intersection(), 1);
        }

        // Periodically store SPPM image in film and write image:
        uint64_t Np = (uint64_t)(iter + 1) * (uint64_t)photonsPerIteration;
        for(int i = 0; i < nPixels; ++i)
        {
            const SPPMPixel &pixel = pixels[i];
            Color3f L = pixel.Ld / (iter + 1.f);
            L += pixel.tau / (Np * Pi * pixel.radius * pixel.radius);
            // Color3f L = pixel.tau / (Np * Pi * pixel.radius * pixel.radius);
            // Color3f L = pixel.tau;
            film->SetPixelColor(pixel.pixel, glm::clamp(L, 0.f, 1.f));
        }
    }


    /*
    // std::vector<Point2i> tilePixels = bounds.GetPoints();
    int pIdx = 0;
    for(Point2i pixel : tilePixels)
    {
        SPPMPixel tempPixel = pixels[pIdx];
        // film->SetPixelColor(pixel, Color3f(1.f, 0.f, 0.f));
        film->SetPixelColor(pixel, glm::clamp(tempPixel.Ld, 0.f, 1.f));
        ++pIdx;
    }
    // Draw points in the grid for basic validation:*/

}
