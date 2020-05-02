#include "thinlenscamera.h"
#include "warpfunctions.h"

ThinLensCamera::ThinLensCamera()
    :Camera(), mFocalLength(29.5f), mLensRadius(1.5f), mSampler(Sampler(100, 0))
{}

ThinLensCamera::ThinLensCamera(unsigned int w, unsigned int h)
    :Camera(w, h), mFocalLength(29.5f), mLensRadius(1.5f), mSampler(Sampler(100, 0))
{}

ThinLensCamera::ThinLensCamera(unsigned int w, unsigned int h, const Vector3f &e, const Vector3f &r, const Vector3f &worldUp)
    :Camera(w, h, e, r, worldUp), mFocalLength(29.5f), mLensRadius(1.5f), mSampler(Sampler(100, 0))
{
    this->RecomputeAttributes();
}

ThinLensCamera::ThinLensCamera(const ThinLensCamera &c): mSampler(Sampler(100, 0))
{
    fovy = c.fovy;
    width = c.width;
    height = c.height;
    near_clip = c.near_clip;
    far_clip = c.far_clip;
    aspect = c.aspect;
    eye = c.eye;
    ref = c.ref;
    look = c.look;
    up = c.up;
    right = c.right;
    world_up = c.world_up;
    V = c.V;
    H = c.H;
    mFocalLength = 29.5f;
    mLensRadius = 1.5f;
    this->RecomputeAttributes();
}

void ThinLensCamera::RecomputeAttributes()
{
    look = glm::normalize(ref - eye);
    right = glm::normalize(glm::cross(look, world_up));
    up = glm::cross(right, look);

    float tan_fovy = tan(glm::radians(fovy/2));
    float len = this->mFocalLength;
    aspect = width/(float)height;
    V = up*len*tan_fovy;
    H = right*len*aspect*tan_fovy;
}

Ray ThinLensCamera::Raycast(const Point2f &pt)
{
    return Raycast(pt.x, pt.y);
}

Ray ThinLensCamera::Raycast(float x, float y)
{
    float ndc_x = (2.f*x/width - 1);
    float ndc_y = (1 - 2.f*y/height);
    return RaycastNDC(ndc_x, ndc_y);
}

Ray ThinLensCamera::RaycastNDC(float ndc_x, float ndc_y)
{
    Point3f focalRef = eye + this->look * this->mFocalLength;
    glm::vec3 pFocal = focalRef + ndc_x*H + ndc_y*V;
    Point2f pLenSample = mSampler.Get2D();
    Point3f pLens = mLensRadius * WarpFunctions::squareToDiskConcentric(pLenSample);
    pLens = eye + up * pLens.y + right * pLens.x;
    Ray result(pLens, glm::normalize(pFocal - pLens));
    // Ray result(eye, glm::normalize(P - eye));

    return result;
}

void ThinLensCamera::create()
{
    std::vector<glm::vec3> pos;
    std::vector<glm::vec3> col;
    std::vector<GLuint> idx;

    //0: Eye position
    pos.push_back(eye);
    //1 - 4: Near clip
        //Lower-left
        Ray r = this->RaycastNDC(-1,-1);
        pos.push_back(eye + r.direction * near_clip);
        //Lower-right
        r = this->RaycastNDC(1,-1);
        pos.push_back(eye + r.direction * near_clip);
        //Upper-right
        r = this->RaycastNDC(1,1);
        pos.push_back(eye + r.direction * near_clip);
        //Upper-left
        r = this->RaycastNDC(-1,1);
        pos.push_back(eye + r.direction * near_clip);
    //5 - 8: Far clip
        //Lower-left
        r = this->RaycastNDC(-1,-1);
        pos.push_back(eye + r.direction * far_clip);
        //Lower-right
        r = this->RaycastNDC(1,-1);
        pos.push_back(eye + r.direction * far_clip);
        //Upper-right
        r = this->RaycastNDC(1,1);
        pos.push_back(eye + r.direction * far_clip);
        //Upper-left
        r = this->RaycastNDC(-1,1);
        pos.push_back(eye + r.direction * far_clip);

    for(int i = 0; i < 9; i++){
        col.push_back(glm::vec3(1,1,1));
    }

    //Frustum lines
    idx.push_back(1);idx.push_back(5);
    idx.push_back(2);idx.push_back(6);
    idx.push_back(3);idx.push_back(7);
    idx.push_back(4);idx.push_back(8);
    //Near clip
    idx.push_back(1);idx.push_back(2);
    idx.push_back(2);idx.push_back(3);
    idx.push_back(3);idx.push_back(4);
    idx.push_back(4);idx.push_back(1);
    //Far clip
    idx.push_back(5);idx.push_back(6);
    idx.push_back(6);idx.push_back(7);
    idx.push_back(7);idx.push_back(8);
    idx.push_back(8);idx.push_back(5);

    //Camera axis
    pos.push_back(eye); col.push_back(glm::vec3(0,0,1)); idx.push_back(9);
    pos.push_back(eye + look); col.push_back(glm::vec3(0,0,1));idx.push_back(10);
    pos.push_back(eye); col.push_back(glm::vec3(1,0,0));idx.push_back(11);
    pos.push_back(eye + right); col.push_back(glm::vec3(1,0,0));idx.push_back(12);
    pos.push_back(eye); col.push_back(glm::vec3(0,1,0));idx.push_back(13);
    pos.push_back(eye + up); col.push_back(glm::vec3(0,1,0));idx.push_back(14);

    count = idx.size();

    bufIdx.create();
    bufIdx.bind();
    bufIdx.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufIdx.allocate(idx.data(), count * sizeof(GLuint));

    bufPos.create();
    bufPos.bind();
    bufPos.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufPos.allocate(pos.data(), pos.size() * sizeof(glm::vec3));

    bufCol.create();
    bufCol.bind();
    bufCol.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufCol.allocate(col.data(), col.size() * sizeof(glm::vec3));
}
