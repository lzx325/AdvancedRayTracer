#ifndef SURFACE_H
#define SURFACE_H
#include <armadillo>
#include <cassert>
#include <memory>
#include "utils.h"
struct LightBase{
    arma::fvec3 position;
    arma::fvec3 intensity;
    LightBase(const arma::fvec3& pos, const arma::fvec3& inten):position(pos),intensity(inten){}
    virtual ~LightBase() {}
};
struct PointLight:public LightBase {
    PointLight(const arma::fvec3& pos, const arma::fvec3& inten):LightBase{pos,inten}{}
    void transform(arma::fmat44 M){
        position=(M*position);
    }
};
struct RectangleLight:public LightBase{
    RectangleLight(const arma::fvec3& pos, const arma::fvec3& inten,const arma::fvec3& cb, const arma::fvec3& ca):LightBase{pos,inten}, cb(cb), ca(ca){}
    arma::fvec3 cb;
    arma::fvec3 ca;
};

struct DRParams{
    arma::fvec2 rectlight_params;
    float motion_t;
};

struct Material {
    Material(float a, float d,const arma::fvec3& dc, float s, float p):ka(a),kd(d),ks(s),p(p),diffuse_color(dc),\
    km(-1),kr(-1),nt(-1),attenuation(arma::fill::zeros){}

    Material(float a, float d,const arma::fvec3& dc, float s, float p,float m):ka(a),kd(d),ks(s),p(p),diffuse_color(dc),km(m),\
    kr(-1),nt(-1),attenuation(arma::fill::zeros){}

    Material(float a, float r, float nt, const arma::fvec3& att):ka(a),kr(r),nt(nt),attenuation(att),\
    kd(-1),ks(-1),km(-1),p(-1),diffuse_color(arma::fill::zeros){}



    Material(): p(0){}
    float ka;
    float kd;
    float ks;
    float km;
    float kr;
    float nt;
    float p;
    arma::fvec3 diffuse_color;
    arma::fvec3 attenuation;
};
struct Ray{
    Ray(const arma::fvec3& ori,const arma::fvec3& dir):origin(ori),direction(dir){}
    Ray()=default;
    arma::fvec3 origin;
    arma::fvec3 direction;
};
struct HitRecord{
    HitRecord(const arma::fvec3& p_, float t_, const arma::fvec3& N_, const Material& mate_,bool out):p(p_),t(t_),N(N_),material(mate_),outside(out){}
    HitRecord():t(std::numeric_limits<float>::max()){}
    arma::fvec3 p;
    float t;
    arma::fvec3 N;
    Material material;
    bool outside;
};

class Surface{
public:
    virtual bool hit(const Ray& ray, float t0, float t1, HitRecord& rec) const =0 ;
    virtual void transform(const arma::fmat44& m) =0;
    virtual Surface* clone()=0;
    Surface()=default;
    Surface(std::string desc):description(desc){}
    Surface(const Surface& rhs)=default;
    Surface(Surface&& rhs)=default;
    Surface& operator=(const Surface& rhs)=default;
    Surface& operator=(Surface&& rhs)=default;
    virtual ~Surface()=default;
    std::string description;
};

struct Sphere: public Surface{
    Sphere(const arma::fvec3 &c, const float r, const Material &m, std::string desc) :Surface(desc), center(c), radius(r), material(m),motion_blur(false) {}
    Sphere(const arma::fvec3 &c, const float r, const Material &m, const arma::fvec3& velocity, std::string desc) :Surface(desc), center(c), radius(r), material(m),motion_blur(true),velocity(velocity) {}
    virtual bool hit(const Ray& ray, float t0, float t1, HitRecord& rec) const override{
        const arma::fvec3& d=ray.direction;
        const arma::fvec3& e=ray.origin;
        const arma::fvec3& c=center;
        assert(fabsf(arma::norm(d)-1.f)<1e-4f);
        assert(t0>=0 && t1>=t0);
        float first=arma::dot(d,e-c);
        float second=first*first;
        float third=arma::dot(e-c,e-c);
        float discrim=second-third+radius*radius;
        if (discrim>=0) {
            float t=-first-sqrtf(discrim);
            if (t<t0) t=-first+sqrtf(discrim);
            if (t<t0 || t>t1) return false;
            arma::fvec3 p=e+t*d;
            arma::fvec3 N=(p-c);
            N/=arma::norm(N);
            bool outside=true;
            if(arma::norm(ray.origin-center)<radius){
                N=-N;
                outside=false;
            }
            p+=1e-3f*N; // move towards normal for a bit
            rec=HitRecord(p,t,N,material,outside);
            return true;
        } else{
            return false;
        }

    }
    Sphere* motion(float t){
        assert(motion_blur && t>=0);
        return new Sphere(center+velocity*t,radius,material,description);
    }
    virtual void transform(const arma::fmat44&) override{
        assert(false);
    }
    virtual Sphere* clone() override{
        assert(false);
        return nullptr;
    }
    arma::fvec3 center;
    float radius;
    Material material;
    bool motion_blur;
    arma::fvec3 velocity;
};
struct Ellipsoid: public Surface{
    Ellipsoid(const arma::fvec3&center, const float a,const float b, const float c, const Material &m, std::string desc=""):\
        Surface(desc),
        center(center),Q(arma::diagmat(arma::fvec4{1/(a*a),1/(b*b),1/(c*c),-1})),material(m){
        arma::fmat44 M_inv=arma::inv(translation_matrix(center));
        this->Q=M_inv.t()*this->Q*M_inv;
    }
    virtual bool hit(const Ray& ray, float t0, float t1, HitRecord& rec) const override{
        assert(fabsf(arma::norm(ray.direction)-1.f)<1e-4f);
        assert(t0>=0 && t1>=t0);

        arma::fvec4 p=to_arma_fvec4(ray.origin,1);
        arma::fvec4 u=to_arma_fvec4(ray.direction,0);
        arma::fmat AA=u.t()*Q*u;
        arma::fmat BB=2*u.t()*Q*p;
        arma::fmat CC=p.t()*Q*p;

        float A=AA(0,0);
        float B=BB(0,0);
        float C=CC(0,0);
        float discrim=B*B-4*A*C;
        if(discrim>=0){
            float t=(-B-sqrtf(discrim))/(2*A);
            if (t<t0){
                t=(-B+sqrtf(discrim))/(2*A);
            }
            if (t<t0||t>t1){
                return false;
            }
            arma::fvec3 hit_point=ray.origin+t*ray.direction;
            arma::fvec4 hit_point_fvec4=to_arma_fvec4(hit_point);
            arma::fvec4 normal_fvec4=Q*hit_point_fvec4;
            arma::fvec3 N=to_arma_fvec3(normal_fvec4);
            N/=arma::norm(N);
            arma::fmat isin=p.t()*Q*p;
            bool outside=true;
            if(isin(0)<0){
                N=-N;
                outside=false;
            }
//            if(arma::dot(ray.direction,N)>=0){
//                std::cout<<ray.origin<<std::endl;
//                std::cout<<ray.direction<<std::endl;
//                std::cout<<N<<std::endl;
//                std::cout<<isin<<std::endl;
//                std::cout<<Q<<std::endl;
//                std::cout<<normal_fvec4<<std::endl;
//                std::cout<<rec.material.diffuse_color<<std::endl;

//            }
            hit_point+=1e-3f*N;
            rec=HitRecord(hit_point,t,N,material,outside);
            return true;
        } else{
            return false;
        }
    }
    virtual void transform(const arma::fmat44& M) override {
        center=to_arma_fvec3(M*to_arma_fvec4(center));
        arma::fmat44 M_inv=arma::inv(M);
        Q=M_inv.t()*Q*M_inv;
        std::cout<<Q<<std::endl;
    }
    virtual Ellipsoid* clone() override{
        return new Ellipsoid(*this);
    }
    arma::fvec3 center;
    arma::fmat44 Q;
    Material material;
};

struct Plane: public Surface{
    Plane():unnormalized_normal(arma::fvec3{1.f,0.f,0.f}),D(0.f){}
    Plane(float A,float B,float C,float D, const Material & material=Material(),std::string desc=""):Surface(desc),unnormalized_normal(arma::fvec3{A,B,C}),D(D),material(material){}
    Plane(const arma::fvec3& u_normal,float D, const Material & material=Material(),std::string desc=""):Surface(desc),unnormalized_normal(u_normal),D(D),material(material){}
    Plane(const arma::fvec4& params, const Material & material=Material(),std::string desc=""):Surface(desc),unnormalized_normal(arma::fvec3{params[0],params[1],params[2]}),D(params[3]),material(material){}
    virtual bool hit(const Ray &ray, float t0, float t1, HitRecord &rec) const override{
        const arma::fvec3& un=unnormalized_normal;
        const arma::fvec3& d=ray.direction;
        const arma::fvec3& e=ray.origin;
        assert(fabsf(arma::norm(d)-1.f)<1e-4f);
        assert(t0>=0 && t1>=t0);
        if (fabsf(arma::dot(un,d))>1e-9f){
            float t=(-arma::dot(un,e)-D)/arma::dot(un,d);
            if (t<t0||t>t1) return false;
            arma::fvec3 p=e+t*d;
            arma::fvec3 N=un/arma::norm(un);
            if (arma::dot(un,e)+D<0){
                N=-N;
            }
            p+=1e-3f*N;
            rec=HitRecord(p,t,N,material,true);
            return true;
        } else{
            return false;
        }
    }
    virtual void transform(const arma::fmat44& M) override{
        const arma::fvec3& un=unnormalized_normal;
        arma::fvec4 n_original={un[0],un[1],un[2],D};
        arma::fmat44 M_inv=arma::inv(M);
        arma::fvec4 n_new=M_inv.t()*n_original;
        unnormalized_normal=to_arma_fvec3(n_new);
        D=n_new(3);
        std::cout<<unnormalized_normal<<std::endl;
        std::cout<<D<<std::endl;
    }
    virtual Plane* clone() override{
        return new Plane(*this);
    }
    arma::fvec3 unnormalized_normal;
    float D;
    Material material;
};
struct Group: public Surface{
    Group():Surface(){}
    Group(const std::vector<std::shared_ptr<Surface>> & surface_group):Surface(), surface_group(surface_group){}
    void add(std::shared_ptr<Surface> surf){
        surface_group.push_back(surf);
    }
    Group(const Group& rhs):Surface(rhs){
        for(size_t i=0;i<rhs.surface_group.size();i++){
            auto p=std::shared_ptr<Surface>(rhs.surface_group[i]->clone());
            this->surface_group.push_back(p);
        }
    }
    Group& operator=(Group rhs){
        swap(*this,rhs);
        return *this;
    }
    friend void swap(Group& lhs,Group&rhs);
    virtual bool hit(const Ray &ray, float t0, float t1, HitRecord &rec) const override{
        assert(t0>=0 && t1>=t0);
        HitRecord rec1;
        bool hit=false;
        for(size_t i=0;i<surface_group.size();i++){
            auto surf=surface_group[i];
            HitRecord rec_tmp;
            if(surf->hit(ray,t0,t1,rec_tmp) && rec_tmp.t<rec1.t){
                hit=true;
                rec1=rec_tmp;
            }
        }
        rec=rec1;
        return hit;
    }
    virtual bool hit(const Ray &ray, float t0, float t1, HitRecord &rec,float time) const {
        assert(t0>=0 && t1>=t0 && time>=0);
        std::vector<std::shared_ptr<Surface> > moved_surface_group;
        moved_surface_group.reserve(surface_group.size());
        for(size_t i=0;i<surface_group.size();i++){
            auto current_ptr=surface_group[i];
            if(typeid(*current_ptr)==typeid(Sphere)){
                auto sphere_ptr=(Sphere*) current_ptr.get();
                if(sphere_ptr->motion_blur){
                    moved_surface_group.push_back(std::shared_ptr<Sphere>(sphere_ptr->motion(time)));
                } else{
                    moved_surface_group.push_back(current_ptr);
                }
            } else{
                moved_surface_group.push_back(current_ptr);
            }
        }
        HitRecord rec1;
        bool hit=false;
        for(size_t i=0;i<moved_surface_group.size();i++){
            auto surf=moved_surface_group[i];
            HitRecord rec_tmp;
            if(surf->hit(ray,t0,t1,rec_tmp) && rec_tmp.t<rec1.t){
                hit=true;
                rec1=rec_tmp;
            }
        }
        rec=rec1;
        return hit;
    }
    virtual int hit_idx(const Ray& ray,float t0,float t1,HitRecord &rec) const {
        assert(t0>=0 && t1>=t0);
        HitRecord rec1;
        int idx=-1;
        for(size_t i=0;i<surface_group.size();i++){
            auto surf=surface_group[i];
            HitRecord rec_tmp;
            if(surf->hit(ray,t0,t1,rec_tmp) && rec_tmp.t<rec1.t){
                idx=int(i);
                rec1=rec_tmp;
            }
        }
        rec=rec1;
        return idx;
    }
    virtual void transform(const arma::fmat44& M) override{
        for(size_t i=0;i<surface_group.size();i++){
            surface_group[i]->transform(M);
        }
    }
    virtual Group* clone() override{
        return new Group(*this);
    }

    std::vector<std::shared_ptr<Surface> > surface_group;
};

inline void swap(Surface&lhs,Surface&rhs){
    using std::swap;
    swap(lhs.description,rhs.description);
}
inline void swap(Group& lhs, Group& rhs){
    using std::swap;
    swap(static_cast<Surface&>(lhs), static_cast<Surface&>(rhs));
    swap(lhs.surface_group,rhs.surface_group);
}
#endif // SURFACE_H
