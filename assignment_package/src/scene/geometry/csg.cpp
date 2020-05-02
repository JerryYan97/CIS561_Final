#include "csg.h"

Intersection CSGNode::Sample(const Point2f &xi, Float *pdf) const
{
    return Intersection();
}

void CSGNode::ComputeTBN(const Point3f &P, Normal3f *nor, Vector3f *tan, Vector3f *bit) const
{
    *nor = glm::normalize(transform.invTransT() * Normal3f(0,0,1));
    CoordinateSystem(*nor, tan, bit);
}

Point2f CSGNode::GetUVCoordinates(const Point3f &point) const
{
    return Point2f(point.x + 0.5f, point.y + 0.5f);
}

float CSGNode::Area() const
{
    return 1.f;
}

void CSGNode::create()
{}

bool segContinuous(CSGSegment& s1, CSGSegment& s2)
{
    if((s1.rightT < s2.leftT) || (s2.rightT < s1.leftT))
    {
        return false;
    }
    else
    {
        return true;
    }
}

CSGSegment CSGSegMerge(CSGSegment& s1, CSGSegment& s2, CSGNodeType Opt)
{
    if(Opt == UNION)
    {
        CSGSegment ans;

        // Get smallest intersection:
        if(s1.leftT < s2.leftT)
        {
            ans.leftT = s1.leftT;
            ans.leftInsect = s1.leftInsect;
        }
        else
        {
            ans.leftT = s2.leftT;
            ans.leftInsect = s2.leftInsect;
        }

        // Get largest intersection:
        if(s1.rightT > s2.rightT)
        {
            ans.rightT = s1.rightT;
            ans.rightInsect = s1.rightInsect;
        }
        else
        {
            ans.rightT = s2.rightT;
            ans.rightInsect = s2.rightInsect;
        }

        return ans;
    }
    else if(Opt == INTERSECT)
    {
        CSGSegment ans;

        // Get largest leftT:
        if(s1.leftT < s2.leftT)
        {
            ans.leftT = s2.leftT;
            ans.leftInsect = s2.leftInsect;
        }
        else
        {
            ans.leftT = s1.leftT;
            ans.leftInsect = s1.leftInsect;
        }

        // Get smallest rightT:
        if(s1.rightT > s2.rightT)
        {
            ans.rightT = s2.rightT;
            ans.rightInsect = s2.rightInsect;
        }
        else
        {
            ans.rightT = s1.rightT;
            ans.rightInsect = s1.rightInsect;
        }

        return ans;
    }
    else if(Opt == DIFFER)
    {
        // Cut s2 from s1:
        if(s1.leftT == s1.rightT)
        {
            // s1 has been totally cut.
            CSGSegment ans;
            ans = s1;
            return ans;
        }

        // s1 still has residual and s1, s2 are continuous:
        if((s1.leftT > s2.leftT) && (s1.rightT < s2.rightT))
        {
            // s2 totally includes s1:
            CSGSegment ans;
            ans.leftT = s1.leftT;
            ans.rightT = s1.leftT;
            return ans;
        }
        else if((s2.rightT > s1.leftT) && (s2.rightT < s1.rightT))
        {
            // s2 intersect s1 from s1's left:
            CSGSegment ans;
            ans.leftT = s2.rightT;
            ans.leftInsect = s2.rightInsect;
            return ans;
        }
        else if((s2.leftT > s1.leftT) && (s2.leftT < s1.rightT))
        {
            // s2 intersect s1 from s1's right:
            CSGSegment ans;
            ans.leftT = s1.leftT;
            ans.leftInsect = s1.leftInsect;
            ans.rightT = s2.leftT;
            ans.rightInsect = s2.leftInsect;
            return ans;
        }
    }

    return CSGSegment();
}

void InputSegment(std::list<CSGSegment>& resSegList, CSGSegment& s, CSGNodeType Opt)
{
    if(resSegList.size() == 0)
    {
        resSegList.push_back(s);
    }
    else
    {
        if(!segContinuous(resSegList.back(), s))
        {
            // Not continuous:
            resSegList.push_back(s);
        }
        else
        {
            // Continuous:
            CSGSegment mergedSeg = CSGSegMerge(resSegList.back(), s, Opt);
            resSegList.pop_back();
            resSegList.push_back(mergedSeg);
        }
    }
}

void UnionOperator(std::list<CSGSegment>& LL, std::list<CSGSegment>& RL, std::list<CSGSegment>& iList)
{
    // Sort segments according to their first intersection's distance:
    if(LL.empty() && RL.empty())
    {
        return;
    }
    else if(LL.empty() && !RL.empty())
    {
        iList = RL;
        return;
    }
    else if(!LL.empty() && RL.empty())
    {
        iList = LL;
        return;
    }

    auto itrL = LL.begin();
    auto itrR = RL.begin();

    while(itrL != LL.end() && itrR != RL.end())
    {
        if(itrL == LL.end())
        {
            for(;itrR != RL.end();itrR++)
            {
                InputSegment(iList, *itrR, UNION);
            }
        }

        if(itrR == RL.end())
        {
            for(;itrL != LL.end();itrL++)
            {
                InputSegment(iList, *itrL, UNION);
            }
        }

        if(itrL->leftT < itrR->leftT)
        {
            InputSegment(iList, *itrL, UNION);
            itrL++;
        }
        else
        {
            InputSegment(iList, *itrR, UNION);
            itrR++;
        }
    }
}

void IntersectOperator(std::list<CSGSegment>& LL, std::list<CSGSegment>& RL, std::list<CSGSegment>& iList)
{
    if(LL.empty() || RL.empty())
    {
        return;
    }
    for(auto itrL = LL.begin(); itrL != LL.end(); itrL++)
    {
        for(auto itrR = RL.begin(); itrR != RL.end(); itrR++)
        {
            if(segContinuous(*itrL, *itrR))
            {
                CSGSegment tempSeg = CSGSegMerge(*itrL, *itrR, INTERSECT);
                iList.push_back(tempSeg);
            }
        }
    }
}

// Left tree segments deduct right tree segments:
void DifferOperator(std::list<CSGSegment>& LL, std::list<CSGSegment>& RL, std::list<CSGSegment>& iList)
{
    if(LL.empty())
    {
        return;
    }
    if(RL.empty())
    {
        iList = LL;
        return;
    }

    // return;
    // Both LL and RL are not empty.
    for(auto itrL = LL.begin(); itrL != LL.end(); itrL++)
    {
        CSGSegment tempSeg = *itrL;
        for(auto itrR = RL.begin(); itrR != RL.end(); itrR++)
        {
            if((itrR->leftT > tempSeg.rightT) || (tempSeg.leftT == tempSeg.rightT))
            {
                break;
            }
            if(segContinuous(tempSeg, *itrR))
            {
                if((tempSeg.leftT < itrR->leftT) && (tempSeg.rightT > itrR->rightT))
                {
                    // itrR is totally included in tempSeg:
                    // Produce two segments:
                    CSGSegment rightTempSeg;
                    rightTempSeg.leftT = itrR->rightT;
                    rightTempSeg.leftInsect = itrR->rightInsect;
                    rightTempSeg.rightT = tempSeg.rightT;
                    rightTempSeg.rightInsect = tempSeg.rightInsect;
                    itrL++;
                    LL.insert(itrL, rightTempSeg);
                    itrL--;
                    itrL--;
                    tempSeg.rightT = itrR->leftT;
                    tempSeg.rightInsect = itrR->leftInsect;
                }
                else
                {
                    tempSeg = CSGSegMerge(tempSeg, *itrR, DIFFER);
                }
            }
        }
        if(tempSeg.leftT != tempSeg.rightT)
        {
            iList.push_back(tempSeg);
        }
    }
}

void CombineSegments(std::list<CSGSegment>& LL, std::list<CSGSegment>& RL, std::list<CSGSegment>& iList, CSGNodeType Opt)
{
    if(Opt == UNION)
    {
        UnionOperator(LL, RL, iList);
    }
    else if(Opt == INTERSECT)
    {
        IntersectOperator(LL, RL, iList);
    }
    else if(Opt == DIFFER)
    {
        DifferOperator(LL, RL, iList);
    }
}

void CSGNode::rayCSGIntersect(const Ray &ray, std::list<CSGSegment>& iList)
{
    // If current node is a geometry, then we need to find out all the intersections between this node
    // and the input ray.
    if(mType == GEOMETRY)
    {
        Intersection tempIntersection;

        if(!selfGeo->Intersect(ray, &tempIntersection))
        {
            return;
        }

        // Discard the intersection that cannot construct a segment:
        // This intersection must have two positive t;
        if(tempIntersection.CSGT1 < 0 || tempIntersection.CSGT2 < 0)
        {
            return;
        }
        else
        {
            // Construct Segment:
            CSGSegment tempSeg;
            tempSeg.leftInsect = tempIntersection;
            tempSeg.rightInsect = tempIntersection;
            tempSeg.leftT = tempIntersection.CSGT1;
            tempSeg.rightT = tempIntersection.CSGT2;
            iList.push_back(tempSeg);
        }
    }
    else
    {
        // An operator should have two children.
        std::list<CSGSegment> LL;
        this->leftChild->rayCSGIntersect(ray, LL);
        if(LL.empty() && mType != UNION)
        {
            return;
        }
        else
        {
            std::list<CSGSegment> RL;
            this->rightChild->rayCSGIntersect(ray, RL);
            // Combine Segments:
            CombineSegments(LL, RL, iList, mType);
            return;
        }
    }
}

// We assume the full binary tree for a CSG tree:
// A node can only have 2 children or no child.
bool CSGNode::Intersect(const Ray &ray, Intersection *isect) const
{
    // return false;
    std::list<CSGSegment> LL;
    this->leftChild->rayCSGIntersect(ray, LL);
    int tempSize = LL.size();
    if(LL.empty() && mType != UNION)
    {
        return false;
    }
    else
    {
        std::list<CSGSegment> RL;
        this->rightChild->rayCSGIntersect(ray, RL);
        // Combine Segments:
        std::list<CSGSegment> tempSegList;
        CombineSegments(LL, RL, tempSegList, mType);
        if(tempSegList.empty())
        {
            return false;
        }
        else
        {
            *isect = tempSegList.begin()->leftInsect;
            return true;
        }
    }
}
