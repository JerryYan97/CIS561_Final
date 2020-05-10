#include "kdtree.h"

KDNode::KDNode()
    : leftChild(nullptr), rightChild(nullptr), axis(0), minCorner(), maxCorner(), particles()
{}

KDNode::~KDNode()
{
    delete leftChild;
    delete rightChild;
}

KDTree::KDTree()
    : root(nullptr)
{}

KDTree::~KDTree()
{
    delete mAllPhotons;
    delete root;
}

// Comparator functions you can use with std::sort to sort vec3s along the cardinal axes
// bool xSort(Photon a, Photon b) { return a.pos.x < b.pos.x; }
// bool ySort(Photon a, Photon b) { return a.pos.y < b.pos.y; }
// bool zSort(Photon a, Photon b) { return a.pos.z < b.pos.z; }
bool xSort(const Photon* a, const Photon* b) { return a->pos.x < b->pos.x; }
bool ySort(const Photon* a, const Photon* b) { return a->pos.y < b->pos.y; }
bool zSort(const Photon* a, const Photon* b) { return a->pos.z < b->pos.z; }

bool xTrisSort(Triangle* a, Triangle* b) { return a->middlePos.x < b->middlePos.x; }
bool yTrisSort(Triangle* a, Triangle* b) { return a->middlePos.y < b->middlePos.y; }
bool zTrisSort(Triangle* a, Triangle* b) { return a->middlePos.z < b->middlePos.z; }


void findCorner(const std::vector<const Photon*> &points, glm::vec3& minCorner, glm::vec3& maxCorner)
{
    minCorner = glm::vec3(points[0]->pos.x, points[0]->pos.y, points[0]->pos.z);
    maxCorner = glm::vec3(points[0]->pos.x, points[0]->pos.y, points[0]->pos.z);
    for(auto itr : points)
    {
        const Photon* tempPtr = itr;
        float x = tempPtr->pos.x;
        float y = tempPtr->pos.y;
        float z = tempPtr->pos.z;
        if(x > maxCorner.x)
        {
            maxCorner.x = x;
        }
        if(x < minCorner.x)
        {
            minCorner.x = x;
        }
        if(y > maxCorner.y)
        {
            maxCorner.y = y;
        }
        if(y < minCorner.y)
        {
            minCorner.y = y;
        }
        if(z > maxCorner.z)
        {
            maxCorner.z = z;
        }
        if(z < minCorner.z)
        {
            minCorner.z = z;
        }
    }
}


int findLeaves(KDNode *v)
{
    int counter = 0;
    if(v->particles.size() == 1)
    {
        return 1;
    }
    else
    {
        if(v->leftChild)
        {
            counter += findLeaves(v->leftChild);
        }
        if(v->rightChild)
        {
            counter += findLeaves(v->rightChild);
        }
        return counter;
    }
}


KDNode* KDTree::buildKdTree(const std::vector<const Photon *> &points,int depth)
{
    if(points.size() == 1)
    {
        KDNode* currNode = new KDNode();
        currNode->minCorner = points[0]->pos;
        currNode->maxCorner = points[0]->pos;
        currNode->particles.push_back(points[0]);
        int depthFlag = depth % 3;
        currNode->axis = depthFlag;
        return currNode;
    }
    else
    {
        std::vector<const Photon *> sortedPoints = points;
        std::vector<const Photon *> setA;
        std::vector<const Photon *> setB;
        glm::vec3 currMinCorner(0.f);
        glm::vec3 currMaxCorner(0.f);
        int depthFlag = depth % 3;
        if(depthFlag == 1)
        {
            // Along x axis.
            std::sort(sortedPoints.begin(), sortedPoints.end(), xSort);
            float xMedian = 0;

            // Get the median along x axis.
            int pSize = sortedPoints.size();
            int sizeFlag = pSize % 2;
            if(sizeFlag == 1)
            {
                // Odd:
                int idx = pSize / 2;
                xMedian = sortedPoints[idx]->pos.x;
            }
            else
            {
                // Even:
                int idx2 = pSize / 2;
                int idx1 = idx2 - 1;
                xMedian = (sortedPoints[idx1]->pos.x + sortedPoints[idx2]->pos.x) / 2.f;
            }

            // Divide points to two sets.
            for(unsigned int i = 0; i < sortedPoints.size(); ++i)
            {
                float currValX = sortedPoints[i]->pos.x;
                if(currValX < xMedian)
                {
                    setA.push_back(sortedPoints[i]);
                }
                else
                {
                    setB.push_back(sortedPoints[i]);
                }
            }
            currMaxCorner = sortedPoints[sortedPoints.size() - 1]->pos;
            currMinCorner = sortedPoints[0]->pos;
        }
        else if(depthFlag == 2)
        {
            // Along y axis.
            std::sort(sortedPoints.begin(), sortedPoints.end(), ySort);
            float yMedian = 0;

            // Get the median along x axis.
            int pSize = sortedPoints.size();
            int sizeFlag = pSize % 2;
            if(sizeFlag == 1)
            {
                // Odd:
                int idx = pSize / 2;
                yMedian = sortedPoints[idx]->pos.y;
            }
            else
            {
                // Even:
                int idx2 = pSize / 2;
                int idx1 = idx2 - 1;
                yMedian = (sortedPoints[idx1]->pos.y + sortedPoints[idx2]->pos.y) / 2.f;
            }

            // Divide points to two sets.
            for(unsigned int i = 0; i < sortedPoints.size(); ++i)
            {
                float currValY = sortedPoints[i]->pos.y;
                if(currValY < yMedian)
                {
                    setA.push_back(sortedPoints[i]);
                }
                else
                {
                    setB.push_back(sortedPoints[i]);
                }
            }
            currMaxCorner = sortedPoints[sortedPoints.size() - 1]->pos;
            currMinCorner = sortedPoints[0]->pos;
        }
        else
        {
            // Along z axis.
            std::sort(sortedPoints.begin(), sortedPoints.end(), zSort);
            float zMedian = 0;

            // Get the median along x axis.
            int pSize = sortedPoints.size();
            int sizeFlag = pSize % 2;
            if(sizeFlag == 1)
            {
                // Odd:
                int idx = pSize / 2;
                zMedian = sortedPoints[idx]->pos.z;
            }
            else
            {
                // Even:
                int idx2 = pSize / 2;
                int idx1 = idx2 - 1;
                zMedian = (sortedPoints[idx1]->pos.z + sortedPoints[idx2]->pos.z) / 2.f;
            }

            // Divide points to two sets.
            for(unsigned int i = 0; i < sortedPoints.size(); ++i)
            {
                float currValZ = sortedPoints[i]->pos.z;
                if(currValZ < zMedian)
                {
                    setA.push_back(sortedPoints[i]);
                }
                else
                {
                    setB.push_back(sortedPoints[i]);
                }
            }
            currMaxCorner = sortedPoints[sortedPoints.size() - 1]->pos;
            currMinCorner = sortedPoints[0]->pos;
        }

        // parent->leftChild = new KDNode();
        // parent->rightChild = new KDNode();
        // depthFlag: 1 -- X, 2 -- Y, 0 -- Z.
        KDNode* currNode = new KDNode();
        // currNode->minCorner = this->minCorner;
        // currNode->maxCorner = this->maxCorner;
        // currNode->minCorner = currMinCorner;
        // currNode->maxCorner = currMaxCorner;
        findCorner(points, currNode->minCorner, currNode->maxCorner);
        currNode->axis = (depthFlag + 2) % 3;
        if(!setA.empty())
        {
            currNode->leftChild = buildKdTree(setA, depth + 1);
        }
        if(!setB.empty())
        {
            currNode->rightChild = buildKdTree(setB, depth + 1);
        }
        return currNode;
    }
}

void KDTree::build(const std::vector<Photon> *points)
{
    mAllPhotons = points;
    if(points->empty())
    {
        return;
    }

    // Build pointers vector:
    std::vector<const Photon*> photonPtrVec(points->size());
    const std::vector<Photon>& photonVec = *mAllPhotons;
    for(unsigned int i = 0; i < points->size(); ++i)
    {
        photonPtrVec[i] = &photonVec[i];
    }

    findCorner(photonPtrVec, minCorner, maxCorner);

    this->root = buildKdTree(photonPtrVec, 1);

    std::cout << "KDTree Build Finished!" << std::endl;
}

void KDTree::buildTris(const QList<Triangle *> &tris)
{
    minCorner = getTriMinCorner(tris);
    maxCorner = getTriMaxCorner(tris);
    root = buildTrisHelper(tris, 0);
    std::cout << "KDTree Build Finished!" << std::endl;
}

glm::vec3 KDTree::getTriMaxCorner(const QList<Triangle *> &tris)
{
    std::vector<float> xarray;
    std::vector<float> yarray;
    std::vector<float> zarray;
    for (auto t : tris)
    {
        xarray.push_back(t->points[0].x);
        xarray.push_back(t->points[1].x);
        xarray.push_back(t->points[2].x);
        yarray.push_back(t->points[0].y);
        yarray.push_back(t->points[1].y);
        yarray.push_back(t->points[2].y);
        zarray.push_back(t->points[0].z);
        zarray.push_back(t->points[1].z);
        zarray.push_back(t->points[2].z);
    }
    auto xresult = std::max_element(xarray.begin(), xarray.end());
    auto yresult = std::max_element(yarray.begin(), yarray.end());
    auto zresult = std::max_element(zarray.begin(), zarray.end());
    return glm::vec3(*xresult, *yresult, *zresult);
}

glm::vec3 KDTree::getTriMinCorner(const QList<Triangle *> &tris)
{
    std::vector<float> xarray;
    std::vector<float> yarray;
    std::vector<float> zarray;
    for (auto t : tris)
    {
        xarray.push_back(t->points[0].x);
        xarray.push_back(t->points[1].x);
        xarray.push_back(t->points[2].x);
        yarray.push_back(t->points[0].y);
        yarray.push_back(t->points[1].y);
        yarray.push_back(t->points[2].y);
        zarray.push_back(t->points[0].z);
        zarray.push_back(t->points[1].z);
        zarray.push_back(t->points[2].z);
    }
    auto xresult = std::min_element(xarray.begin(), xarray.end());
    auto yresult = std::min_element(yarray.begin(), yarray.end());
    auto zresult = std::min_element(zarray.begin(), zarray.end());
    return glm::vec3(*xresult, *yresult, *zresult);
}

glm::vec3 KDTree::getMaxCorner(const std::vector<Photon> &points)
{
    std::vector<float> xarray;
    std::vector<float> yarray;
    std::vector<float> zarray;
    for (unsigned int i = 0; i < points.size(); i++)
    {
        xarray.push_back(points[i].pos.x);
        yarray.push_back(points[i].pos.y);
        zarray.push_back(points[i].pos.z);
    }
    auto xresult = std::max_element(xarray.begin(), xarray.end());
    auto yresult = std::max_element(yarray.begin(), yarray.end());
    auto zresult = std::max_element(zarray.begin(), zarray.end());
    return glm::vec3(*xresult, *yresult, *zresult);
}

glm::vec3 KDTree::getMinCorner(const std::vector<Photon> &points)
{
    std::vector<float> xarray;
    std::vector<float> yarray;
    std::vector<float> zarray;
    for (unsigned int i = 0; i < points.size(); i++)
    {
        xarray.push_back(points[i].pos.x);
        yarray.push_back(points[i].pos.y);
        zarray.push_back(points[i].pos.z);
    }
    auto xresult = std::min_element(xarray.begin(), xarray.end());
    auto yresult = std::min_element(yarray.begin(), yarray.end());
    auto zresult = std::min_element(zarray.begin(), zarray.end());
    return glm::vec3(*xresult, *yresult, *zresult);
}


KDNode* KDTree::buildTrisHelper(const QList<Triangle *> &tris, int depth)
{
    KDNode *mynode = new KDNode();
    QList<Triangle*> sortTris = tris;
    unsigned int axis = depth & 3;
    if (sortTris.size() > 1)
    {
        if (axis == 0)
        {
            std::sort(sortTris.begin(), sortTris.end(), xTrisSort);
        }
        else if (axis == 1)
        {
            std::sort(sortTris.begin(), sortTris.end(), yTrisSort);
        }
        else
        {
            std::sort(sortTris.begin(), sortTris.end(), zTrisSort);
        }
        int median = sortTris.size() / 2;
        QList<Triangle*> leftPoints(sortTris.mid(0, median));
        QList<Triangle*> rightPoints(sortTris.mid(median));
        mynode->axis = axis;
        mynode->minCorner = getTriMinCorner(sortTris);
        mynode->maxCorner = getTriMaxCorner(sortTris);
        mynode->leftChild = buildTrisHelper(leftPoints, depth+1);
        mynode->rightChild = buildTrisHelper(rightPoints, depth+1);
    }
    else
    {
        mynode->axis = axis;
        mynode->triangles = sortTris;
        mynode->minCorner = getTriMinCorner(sortTris);
        mynode->maxCorner = getTriMaxCorner(sortTris);
    }
    return mynode;
}


inline bool inSphere(const glm::vec3& c, float r, const glm::vec3& v)
{
    float rVSquare = (v.x - c.x) * (v.x - c.x) + (v.y - c.y) * (v.y - c.y) + (v.z - c.z) * (v.z - c.z);
    r = r * r;
    if(rVSquare <= r)
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline float squared(const float& v) { return v * v; }
bool doesCubeIntersectSphere(const glm::vec3& C1, const glm::vec3& C2, const glm::vec3& S, const float& R)
{
    float dist_squared = R * R;
    /* assume C1 and C2 are element-wise sorted, if not, do that now */
    if (S.x < C1.x) dist_squared -= squared(S.x - C1.x);
    else if (S.x > C2.x) dist_squared -= squared(S.x - C2.x);
    if (S.y < C1.y) dist_squared -= squared(S.y - C1.y);
    else if (S.y > C2.y) dist_squared -= squared(S.y - C2.y);
    if (S.z < C1.z) dist_squared -= squared(S.z - C1.z);
    else if (S.z > C2.z) dist_squared -= squared(S.z - C2.z);
    return dist_squared > 0;
}


// Between Sphere and cuboid.
inline IntersectState checkIntersection(const glm::vec3& c, const float& r, const glm::vec3& minCorner, const glm::vec3& maxCorner)
{
    std::vector<glm::vec3> cuboidCorners(8);
    std::vector<float> tempX{minCorner.x, maxCorner.x};
    std::vector<float> tempY{minCorner.y, maxCorner.y};
    std::vector<float> tempZ{minCorner.z, maxCorner.z};

    for(int i = 0; i < 2; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            for(int k = 0; k < 2; ++k)
            {
                cuboidCorners[i * 4 + j * 2 + k] = glm::vec3(tempX[i], tempY[j], tempZ[k]);
            }
        }
    }

    // Check whether 8 points are all in sphere.
    int intersectionCounter = 0;
    for(int i = 0; i < 8; ++i)
    {
        if(inSphere(c, r, cuboidCorners[i]))
        {
            intersectionCounter++;
        }
    }

    if(intersectionCounter == 8)
    {
        return FullyContained;
    }
    else
    {
        if(doesCubeIntersectSphere(minCorner, maxCorner, c, r))
        {
            return PartlyIntersect;
        }
        else
        {
            return NoIntersect;
        }
    }
}


void KDTree::reportSubTree(KDNode *v, std::vector<const Photon*> &reportedPoints) const
{
    if(v->particles.size() == 1)
    {
        reportedPoints.push_back(v->particles[0]);
    }
    else
    {
        if(v->leftChild)
        {
            reportSubTree(v->leftChild, reportedPoints);
        }
        if(v->rightChild)
        {
            reportSubTree(v->rightChild, reportedPoints);
        }
    }
}


void KDTree::searchKdTree(KDNode* v, std::vector<const Photon*>& reportedPoints, const glm::vec3& c, const float& r) const
{
    if (v->particles.size() != 0 && (v->leftChild == nullptr && v->rightChild == nullptr))
    {
        if(glm::distance(v->particles[0]->pos, c) <= r)
        {
            reportedPoints.push_back(v->particles[0]);
        }
    }
    else
    {
        if(v->leftChild)
        {
            if(doesCubeIntersectSphere(v->leftChild->minCorner, v->leftChild->maxCorner, c, r))
            {
                searchKdTree(v->leftChild, reportedPoints, c, r);
            }
        }
        if(v->rightChild)
        {
            if(doesCubeIntersectSphere(v->rightChild->minCorner, v->rightChild->maxCorner, c, r))
            {
                searchKdTree(v->rightChild, reportedPoints, c, r);
            }
        }
    }
}

std::vector<const Photon*> KDTree::particlesInSphere(const glm::vec3& c, const float& r) const
{
    //TODO
    std::vector<const Photon*> reportedPoints;
    if(this->root)
    {
        searchKdTree(this->root, reportedPoints, c, r);
    }
    // std::cout << "reportedpoints num: " << reportedPoints.size() << std::endl;
    return reportedPoints;
}

bool KDTree::findIntersectTris(const Ray &r, QList<Triangle *> &tris)
{
    if (hitBoundBox(minCorner, maxCorner, r)){
        traverse(root, r, tris);
        if (tris.size() != 0) return true;
    }
    return false;
}

void KDTree::traverse(KDNode *node, const Ray &r, QList<Triangle*> &tris)
{
    if (node->triangles.size() != 0 && (node->leftChild == nullptr && node->rightChild == nullptr))
    {
        tris.append(node->triangles[0]);
    }
    else
    {
        if (hitBoundBox(node->leftChild->minCorner, node->rightChild->maxCorner, r))
        {
            traverse(node->leftChild, r, tris);
        }
        if (hitBoundBox(node->rightChild->minCorner, node->rightChild->maxCorner, r))
        {
            traverse(node->rightChild, r, tris);
        }
    }
}

bool KDTree::hitBoundBox(const Point3f &minC, const Point3f &maxC, const Ray &r)
{
    Point3f bounds[2];
    bounds[0] = minC;
    bounds[1] = maxC;
    Vector3f invdir = 1.f / r.direction;
    int sign[3];
    sign[0] = (invdir.x < 0);
    sign[1] = (invdir.y < 0);
    sign[2] = (invdir.z < 0);
    float tmin, tmax, tymin, tymax, tzmin, tzmax;

    tmin = (bounds[sign[0]].x - r.origin.x) * invdir.x;
    tmax = (bounds[1-sign[0]].x - r.origin.x) * invdir.x;
    tymin = (bounds[sign[1]].y - r.origin.y) * invdir.y;
    tymax = (bounds[1-sign[1]].y - r.origin.y) * invdir.y;

    if ((tmin > tymax) || (tymin > tmax)) return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    tzmin = (bounds[sign[2]].z - r.origin.z) * invdir.z;
    tzmax = (bounds[1-sign[2]].z - r.origin.z) * invdir.z;

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;

    return true;

}

bool KDTree::intersect(KDNode *node, glm::vec3 c, float r) const
{
    float sumDist = 0;
    float r2 = r*r;
    for (int i = 0; i < 3; i++)
    {
        if (c[i] < node->minCorner[i]) sumDist += glm::length2(c[i] - node->minCorner[i]);
        else if (c[i] > node->maxCorner[i]) sumDist += glm::length2(c[i] - node->maxCorner[i]);
    }
    return sumDist <= r2;
}

void KDTree::clear()
{
    delete root;
    root = nullptr;
}
