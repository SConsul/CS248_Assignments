
#include "../scene/skeleton.h"
#include<iostream>
Vec3 closest_on_line_segment(Vec3 start, Vec3 end, Vec3 point) {

    // TODO(Animation): Task 3
    // Return the closest point to 'point' on the line segment from start to end
    
    Vec3 s2e = end - start, s2p = point - start;
    float length = s2e.norm();
    Vec3 s2eNorm = s2e.normalize();
    
    float projection = dot(s2p, s2eNorm);
    if(projection >= length){return end;}
    else if(projection <= 0){return start;}
    else{
        return projection * s2eNorm + start;
    }
}

Mat4 Joint::joint_to_bind() const {

    // TODO(Animation): Task 2

    // Return a matrix transforming points in the space of this joint
    // to points in skeleton space in bind position.

    // Bind position implies that all joints have pose = Vec3{0.0f}

    // You will need to traverse the joint heirarchy. This should
    // not take into account Skeleton::base_pos
    // return Mat4::I;
    Mat4 T = Mat4::I;
    for (Joint* j=this->parent; j!=nullptr; j=j->parent){
        T = Mat4::translate(j->extent)*T;
    }
    return T;

}

Mat4 Joint::joint_to_posed() const {

    // TODO(Animation): Task 2

    // Return a matrix transforming points in the space of this joint
    // to points in skeleton space, taking into account joint poses.

    // You will need to traverse the joint heirarchy. This should
    // not take into account Skeleton::base_pos
    Mat4 T = Mat4::euler(this->pose);
    for (Joint* j=this->parent; j!=nullptr; j=j->parent){
        T = Mat4::translate(j->extent)*T;
        T = Mat4::euler(j->pose)*T;
    }
    return T;
}

Vec3 Skeleton::end_of(Joint* j) {

    // TODO(Animation): Task 2

    // Return the bind position of the endpoint of joint j in object space.
    // This should take into account Skeleton::base_pos.
    return j->joint_to_bind()*j->extent + this->base_pos;
}

Vec3 Skeleton::posed_end_of(Joint* j) {

    // TODO(Animation): Task 2

    // Return the posed position of the endpoint of joint j in object space.
    // This should take into account Skeleton::base_pos.
    return j->joint_to_posed()*j->extent + this->base_pos;
}

Mat4 Skeleton::joint_to_bind(const Joint* j) const {

    // TODO(Animation): Task 2

    // Return a matrix transforming points in joint j's space to object space in
    // bind position. This should take into account Skeleton::base_pos.
    Mat4 T = Mat4::I;
    for (Joint* i=j->parent; i!=nullptr; i=i->parent){
        T = Mat4::translate(i->extent)*T;
        T = Mat4::euler(i->pose)*T;
    }

    T = Mat4::translate(this->base_pos) * T;
    return T;
}

Mat4 Skeleton::joint_to_posed(const Joint* j) const {

    // TODO(Animation): Task 2

    // Return a matrix transforming points in joint j's space to object space with
    // poses. This should take into account Skeleton::base_pos.
    Mat4 T = Mat4::I;
    for (Joint* i=j->parent; i!=nullptr; i=i->parent){
        T = Mat4::translate(i->extent)*T;
    }

    T = Mat4::translate(this->base_pos) * T;
    return T;
}

void Skeleton::find_joints(const GL::Mesh& mesh,
                           std::unordered_map<unsigned int, std::vector<Joint*>>& map) {

    // TODO(Animation): Task 3

    // Construct a mapping from vertex indices to lists of joints in this skeleton
    // that should effect the vertex at that index. A joint should effect a vertex
    // if it is within Joint::radius distance of the bone's line segment in bind position.

    const std::vector<GL::Mesh::Vert>& verts = mesh.verts();

    // For each i in [0, verts.size()), map[i] should contain the list of joints that
    // effect vertex i. Note that i is NOT Vert::id! i is the index in verts.

    for_joints([&](Joint* j) {
        // What vertices does joint j effect?
        for(size_t i=0; i< verts.size(); i++){
            Vec3 closestPoint = closest_on_line_segment(base_of(j), end_of(j), verts[i].pos /*+ this->base_pos*/);
            float distance = (closestPoint - verts[i].pos).norm();
            if(distance < j->radius){
                map[i].push_back(j);
            }
        }
    });
}

void Skeleton::skin(const GL::Mesh& input, GL::Mesh& output,
                    const std::unordered_map<unsigned int, std::vector<Joint*>>& map) {

    // TODO(Animation): Task 3

    // Apply bone poses & weights to the vertices of the input (bind position) mesh
    // and store the result in the output mesh. See the task description for details.
    // map was computed by find_joints, hence gives a mapping from vertex index to
    // the list of bones the vertex should be effected by.

    // Currently, this just copies the input to the output without modification.

    std::vector<GL::Mesh::Vert> verts = input.verts();
    for(size_t i = 0; i < verts.size(); i++) {
        // Skin vertex i. Note that its position is given in object bind space.
        if(map.count(i) <= 0){continue;}
        Vec3 num;
        float den = 0.0;
        for(Joint* j : map.at(i)){
            Vec3 closestPoint = closest_on_line_segment(base_of(j), end_of(j), verts[i].pos /*+ this->base_pos*/);
            float dist_ij_inv = 1.0/((closestPoint - (verts[i].pos /*+ this->base_pos*/)).norm());

            Vec3 v_ij = (j->joint_to_posed() * j->joint_to_bind().inverse() * (verts[i].pos - this->base_pos)) + this->base_pos;
            
            num+= v_ij * dist_ij_inv;
            den+= dist_ij_inv;
        }
        verts[i].pos = num/den;
        
    }

    std::vector<GL::Mesh::Index> idxs = input.indices();
    output.recreate(std::move(verts), std::move(idxs));
}

void Joint::compute_gradient(Vec3 target, Vec3 current) {

    // TODO(Animation): Task 2

    // Computes the gradient of IK energy for this joint and, should be called
    // recursively upward in the heirarchy. Each call should storing the result
    // in the angle_gradient for this joint.

    // Target is the position of the IK handle in skeleton space.
    // Current is the end position of the IK'd joint in skeleton space.

    std::vector<Joint*> hier;
    for (Joint* j=this; j!=nullptr; j=j->parent){
        hier.push_back(j);
    }

    Vec3 delta = current-target;
    for(auto iter=hier.rbegin(); iter!=hier.rend(); iter++){
        Joint* j = *iter;
        Mat4 T = j->joint_to_posed();
        Vec3 base = Vec3(); //base of point is origin in that joint space
        Vec3 p = current -  j->joint_to_posed()*base;

        Vec4 X_axis_h = T*Vec4(1,0,0,0);
        Vec3 X_axis = Vec3(X_axis_h.x, X_axis_h.y, X_axis_h.z).normalize();
        j->angle_gradient.x += dot(delta,cross(X_axis,p));

        Vec4 Y_axis_h = T*Vec4(0,1,0,0);
        Vec3 Y_axis = Vec3(Y_axis_h.x, Y_axis_h.y, Y_axis_h.z).normalize();
        j->angle_gradient.y += dot(delta,cross(Y_axis,p));

        Vec4 Z_axis_h = T*Vec4(0,0,1,0);
        Vec3 Z_axis = Vec3(Z_axis_h.x, Z_axis_h.y, Z_axis_h.z).normalize();
        j->angle_gradient.z += dot(delta,cross(Z_axis,p));
    }

}

void Skeleton::step_ik(std::vector<IK_Handle*> active_handles) {
    // TODO(Animation): Task 2
    // Do several iterations of Jacobian Transpose gradient descent for IK
    for(int iter = 0; iter<5  ; iter++){
        //zero_grad()
        for_joints([&](Joint* j){j->angle_gradient = Vec3();});

        //compute gradient
        for(auto handle: active_handles){
            handle->joint->compute_gradient(handle->target, posed_end_of(handle->joint)- this->base_pos);
        }
            
        //step
        float lr = 0.2;
        for_joints([&](Joint* j){j->pose -=lr*j->angle_gradient;});
    }
}
