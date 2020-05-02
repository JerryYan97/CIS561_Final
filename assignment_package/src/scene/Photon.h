#pragma once
#include <la.h>
#include <globals.h>

class Primitive;

class Photon
{
public:
    Point3f pos;
    Color3f color;
    Vector3f wi;
    Vector3f hitNormal;
    Primitive const *hitPrim;
    // QString attachedObjName;
    Photon()
        : pos(Point3f(0.f)), color(Color3f(0.f)), wi(Vector3f(0.f))
    {}
    Photon(Point3f p, Color3f c, Vector3f w)
        : pos(p), color(c), wi(w)
    {}
};
