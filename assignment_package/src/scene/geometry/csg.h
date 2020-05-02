#ifndef CSG_H
#define CSG_H
#include "shape.h"
#include "primitive.h"

enum CSGNodeType{
    UNION, INTERSECT, DIFFER, GEOMETRY
};

class CSGIntersection
{
public:
    CSGIntersection(Intersection& iSect, bool iIn, float iDis)
        :mInsect(iSect), mIn(iIn), mDistance(iDis)
    {}
    Intersection mInsect;
    bool mIn;
    float mDistance;
};

class CSGSegment
{
public:
    // std::vector<CSGIntersection> mInsectVec;
    Intersection leftInsect;
    float leftT;
    Intersection rightInsect;
    float rightT;
};


class CSGNode : public Shape
{
public:
    CSGNode* leftChild;
    CSGNode* rightChild;
    Shape* selfGeo;
    CSGNodeType mType;

    virtual bool Intersect(const Ray &ray, Intersection *isect) const;
    void rayCSGIntersect(const Ray &ray, std::list<CSGSegment>& iList);
    virtual Point2f GetUVCoordinates(const Point3f &point) const;
    virtual void ComputeTBN(const Point3f& P, Normal3f* nor, Vector3f* tan, Vector3f* bit) const;

    virtual float Area() const;

    // Sample a point on the surface of the shape and return the PDF with
    // respect to area on the surface.
    virtual Intersection Sample(const Point2f &xi, Float *pdf) const;

    void create();
};

#endif // CSG_H
