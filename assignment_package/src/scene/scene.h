#pragma once
#include <QList>
#include <raytracing/film.h>
#include "thinlenscamera.h"
// #include <scene/camera.h>
#include <scene/lights/light.h>

class Primitive;
class Material;
class Light;
class Shape;

class Scene
{
public:
    Scene();
    QList<std::shared_ptr<Primitive>> primitives;
    QList<std::shared_ptr<Material>> materials;
    QList<std::shared_ptr<Light>> lights;

    QList<std::shared_ptr<Shape>> CSGShapes;

    Camera* mCameraPtr;
    Camera camera;
    ThinLensCamera tlCamera;


    Film film;

    void SetCamera(const Camera &c);
    void SetThinLensCamera(const ThinLensCamera &c);

    void CreateTestScene();
    void Clear();

    bool Intersect(const Ray& ray, Intersection* isect) const;
};
