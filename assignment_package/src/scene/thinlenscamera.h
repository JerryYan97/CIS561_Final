#ifndef THINLENSCAMERA_H
#define THINLENSCAMERA_H
#include "camera.h"
#include <samplers/sampler.h>


class ThinLensCamera : public Camera
{
public:
    ThinLensCamera();
    ThinLensCamera(unsigned int w, unsigned int h);
    ThinLensCamera(unsigned int w, unsigned int h, const Vector3f &e, const Vector3f &r, const Vector3f &worldUp);
    ThinLensCamera(const ThinLensCamera &c);

    virtual void RecomputeAttributes() override;

    virtual Ray Raycast(const Point2f &pt) override;         //Creates a ray in 3D space given a 2D point on the screen, in screen coordinates.
    virtual Ray Raycast(float x, float y) override;            //Same as above, but takes two floats rather than a vec2.
    virtual Ray RaycastNDC(float ndc_x, float ndc_y) override; //Creates a ray in 3D space given a 2D point in normalized device coordinates.

    void create() override;
    float mFocalLength;
    float mLensRadius;
    Sampler mSampler;
    // Point2f sample;
};

#endif // THINLENSCAMERA_H
