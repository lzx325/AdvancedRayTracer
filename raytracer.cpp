#include "raytracer.h"
#include <cstring>
#include <iostream>
#include <QDebug>
#include <omp.h>
#include <memory>
#include "surface.h"
#include "utils.h"
float R(float nt, float cos_theta){

    assert(cos_theta>=0 && cos_theta<=1);
    float R0=powf((nt-1)/(nt+1),2);
    return R0+(1-R0)*(1-powf(cos_theta,5));
}

arma::fvec3 raycolor(const Ray& ray,const Group& scene, const std::vector<std::shared_ptr<LightBase> >& lights,
               const arma::fvec3& Ia, const arma::fvec3& bg, float t0, float t1, const DRParams& drparams, size_t recursion){
    arma::fvec3 return_color;
    HitRecord rec,srec;
    if(recursion>4 || !scene.hit(ray,t0,t1,rec,drparams.motion_t)){
        return_color=bg;
    } else{
        arma::fvec3 Id{0,0,0},Is{0,0,0},Im{0,0,0},Ir{0,0,0};
        if(rec.material.kd>=0 && rec.material.ks>=0){
            for(size_t i=0;i<lights.size();i++){
                arma::fvec3 position=lights[i]->position;

                if(typeid(*lights[i])==typeid(RectangleLight)){

                    auto rectangle_ptr=(RectangleLight*) lights[i].get();
                    position+=drparams.rectlight_params[0]*rectangle_ptr->ca
                            +drparams.rectlight_params[1]*rectangle_ptr->cb;
                }
                arma::fvec3 l=(position-rec.p);
                l/=arma::norm(l);
                Ray sray(rec.p,l);
                float light_distance=arma::norm(position-rec.p);
                if(!scene.hit(sray,1e-2f,std::numeric_limits<float>::max(),srec,drparams.motion_t) || srec.t>light_distance){
                    arma::fvec3 h=(l-ray.direction);
                    h/=arma::norm(h);
                    Id+=lights[i]->intensity*std::max(0.f,arma::dot(rec.N,l));
                    Is+=lights[i]->intensity*powf(std::max(0.f,arma::dot(rec.N,h)),rec.material.p);
                }
            }
        }

        arma::fvec3 reflect_dir=reflect(ray.direction,rec.N);

        assert(rec.material.km<0||rec.material.kr<0);
        if(rec.material.km>=0){
            Ray reflect_ray(rec.p,reflect_dir);
            Im=raycolor(reflect_ray, scene, lights, Ia, bg, t0, t1,drparams,recursion+1);

        }

        if(rec.material.kr>=0){
            arma::fvec3 refract_dir;

            if(refract(ray.direction,rec.N,rec.outside?rec.material.nt:1/rec.material.nt,refract_dir)){
                // reflection and refraction
                Ray reflect_ray(rec.p,reflect_dir);
                Ray refract_ray(rec.p-2e-3f*rec.N,refract_dir);
                arma::fvec3 reflect_color=raycolor(reflect_ray, scene, lights,
                                                   Ia, bg, t0, t1,drparams,recursion+1);
                arma::fvec3 refract_color=raycolor(refract_ray, scene, lights,
                                                   Ia, bg, t0, t1,drparams,recursion+1);
                float cos_theta=-arma::dot(ray.direction,rec.N);
                if(!(cos_theta>=0 && cos_theta<=1)){
                    std::cout<<cos_theta<<std::endl;
                }
                if(cos_theta<0) cos_theta=0;
                if(cos_theta>1) cos_theta=1;
//                float R_val=R(rec.material.nt,cos_theta);
                float R_val=0.1;
                Ir=R_val*reflect_color+(1-R_val)*refract_color;

            } else{
                Ray reflect_ray(rec.p,reflect_dir);
                Ir=raycolor(reflect_ray, scene, lights, Ia, bg, t0, t1,drparams,recursion+1);
            }
        }
        assert(arma::all(rec.material.attenuation>=0));
        return_color= (rec.material.ka>0?rec.material.ka:0)*Ia%rec.material.diffuse_color
                +(rec.material.kd>0?rec.material.kd:0)*Id%rec.material.diffuse_color
                +(rec.material.ks>0?rec.material.ks:0)*Is
                +(rec.material.km>0?rec.material.km:0)*Im
                +(rec.material.kr>0?rec.material.kr:0)*Ir;

        if(! rec.outside){
            return_color%=arma::exp(-rec.material.attenuation*rec.t);
        }
    }
    return return_color;
}

arma::fvec3 raycolor(const Ray& ray,const Group& scene, const std::vector<std::shared_ptr<LightBase>>& lights,
                     const arma::fvec3& Ia, const arma::fvec3& bg, float t0, float t1, const DRParams& drparams){
    return raycolor(ray, scene, lights, Ia, bg, t0, t1,drparams,0);
}


size_t img_width=1600;
size_t img_height=1200;

std::vector<arma::fvec3> render(const Group & scene, const std::vector<std::shared_ptr<LightBase> > &lights,float tan_vfov=1)
{
    // render options
    size_t antialiasing=5;
    arma::fvec3 viewpoint{0,0,1};
    float l=-tan_vfov*4.f/3.f,r=tan_vfov*4.f/3.f,b=-tan_vfov,t=tan_vfov;
    arma::fvec3 Ia{1.f,1.f,1.f};
    arma::fvec3 bg{0,0,0};
    float t0=0;
    float t1=1000;
    std::vector<arma::fvec3> img(img_width * img_height);

    #pragma omp parallel for
    for (int j = 0; j < (int)img_height; j++)
    {
        arma::fmat r_mat(antialiasing*antialiasing,2);
        arma::fmat s_mat(antialiasing*antialiasing,2);
        arma::fvec m_vec(antialiasing*antialiasing);
        DRParams dummy_drparams{{0,0},0};
        for (size_t i = 0; i < img_width; i++)
        {
            if(antialiasing>1){
                for(size_t p=0;p<antialiasing;p++){
                    for(size_t q=0;q<antialiasing;q++){
                        r_mat(p*antialiasing+q,0)=(p+arma::randu())/antialiasing;
                        r_mat(p*antialiasing+q,1)=(q+arma::randu())/antialiasing;
                        s_mat(p*antialiasing+q,0)=(p+arma::randu())/antialiasing;
                        s_mat(p*antialiasing+q,1)=(q+arma::randu())/antialiasing;
                        m_vec(p*antialiasing+q)=(p*antialiasing+q+arma::randu())/(antialiasing*antialiasing);
                    }
                }
                for(int k=antialiasing*antialiasing-1;k>=1;k--){
                    int l=arma::randi(1,arma::distr_param(0,k))(0);
                    std::swap(s_mat(l,0),s_mat(k,0));
                    std::swap(s_mat(l,1),s_mat(k,1));
                    std::swap(m_vec(l),m_vec(k));
                }
                arma::fvec3 c(arma::fill::zeros);
                for(size_t k=0;k<antialiasing*antialiasing;k++){
                    float dir_x = l+(r-l)*(i+r_mat(k,0))/img_width;
                    float dir_y = t-(t-b)*(j + r_mat(k,1))/img_height; // this flips the image at the same time
                    float dir_z = -1;
                    arma::fvec3 dir={dir_x,dir_y,dir_z};
                    dir/=arma::norm(dir);
                    c+=raycolor(Ray{viewpoint,dir}, scene, lights,
                                Ia, bg, t0, t1, {{s_mat(k,0),s_mat(k,1)},m_vec(k)});
                }

                img[i+j*img_width]=c/(antialiasing*antialiasing);

            } else{
                float dir_x = l+(r-l)*(i+0.5f)/img_width;
                float dir_y = t-(t-b)*(j + 0.5f)/img_height; // this flips the image at the same time
                float dir_z = -1;
                arma::fvec3 dir={dir_x,dir_y,dir_z};
                dir/=arma::norm(dir);
                arma::fvec3 c;
                c=raycolor(Ray{viewpoint,dir}, scene, lights,
                                        Ia, bg, t0, t1,dummy_drparams);
                img[i+j*img_width]=c;
            }


        }
    }
    return img;
}




