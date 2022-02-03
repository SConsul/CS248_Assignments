
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>
#include <iostream>
#include "../geometry/halfedge.h"
#include "debug.h"

Halfedge_Mesh::HalfedgeRef degenerateHE[2];
uint8_t numDegenerateHE = 0;

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it does not want to perform the operation for
    whatever reason (e.g. you don't want to allow the user to erase the last vertex).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementaiton, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/

/*
    This method should replace the given vertex and all its neighboring
    edges and faces with a single face, returning the new face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_vertex(Halfedge_Mesh::VertexRef v) {

    (void)v;
    return std::nullopt;
}

/*
    This method should erase the given edge and return an iterator to the
    merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_edge(Halfedge_Mesh::EdgeRef e) {

    (void)e;
    return std::nullopt;
}

/*
    This method should collapse the given edge and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(Halfedge_Mesh::EdgeRef e) {
    numDegenerateHE = 0;

    // If only 3 vertices are left, do nothing
    if(vertices.size() - verased.size() <=3){
        return std::nullopt;
    }
    
    // std::cout << "Collapsing edge " << e->id() << std::endl;
    //check for number of adjacent vertices
    int num_neigh_v = 0;
    VertexRef v0 = e->halfedge()->vertex(), v1 = e->halfedge()->twin()->vertex();
    HalfedgeRef h0 = e->halfedge()->twin()->next(), h1= e->halfedge()->next();
    while(h0!=e->halfedge()){
        if(!h0->is_boundary() && h0->next()->twin()->vertex()==v1) num_neigh_v++;
        h0 = h0->twin()->next();
    }
    while(h1!=e->halfedge()->twin()){
        if(!h1->is_boundary() && h1->next()->twin()->vertex()==v0) num_neigh_v++;
        h1 = h1->twin()->next();
    }
    if(num_neigh_v>2){
        std::cout<<"edge collapse would lead to non-manifold mesh"<<std::endl;
        return std::nullopt;
    } 
    
    //dangling vertex check
    int len1 = 1, len2=1;
    HalfedgeRef h = e->halfedge();
    while(h->next()!=e->halfedge() && len1<=3){
        h = h->next(); 
        len1++;
    }
    h = e->halfedge()->twin();
    while(h->next()!=e->halfedge()->twin() && len2<=3){
        h = h->next(); 
        len2++;
    }
    if( (len1<=3 && len2<=3) && 
    (e->halfedge()->next()->twin()->vertex()->degree()<=2 ||
     e->halfedge()->twin()->next()->twin()->vertex()->degree()<=2)){
        return std::nullopt;
    }

    struct Changes
    {
        HalfedgeRef h;
        HalfedgeRef h_twin_prev;
    };
    Changes currChanges[2];

    struct ChangesND  // Changes for non-degenerate faces
    {
        HalfedgeRef h, h_prev_prev;
    };
    std::vector<ChangesND> currChangesND;

    HalfedgeRef toBeErased[2];
    
    VertexRef v_new = new_vertex();
    v_new->pos = e->center();
    h = e->halfedge();
    v_new->halfedge() = h->next();

    //iterate over edges of v0 = h->vertex()
    h = e->halfedge()->twin()->next();
    while(h!= e->halfedge()){
        //move HE from v to v_new
        h->set_neighbors(h->next(),h->twin(),v_new, h->edge(), h->face());
        h = h->twin()->next(); 
    }
    // h->set_neighbors(h->next(),h->twin(),v_new, h->edge(), h->face());
    
    //iterate over edges of v1 = h->twin()->vertex()
    h = e->halfedge()->next();
    while(h!= e->halfedge()->twin()){
        //move HE from v to v_new
        h->set_neighbors(h->next(),h->twin(),v_new, h->edge(), h->face());
        h = h->twin()->next(); 
    }
    // h->set_neighbors(h->next(),h->twin(),v_new, h->edge(), h->face());

    // //check for degeneracy while iterating over edeges of v_new
    h = e->halfedge()->next();
    HalfedgeRef comparisonHE = e->halfedge()->next();
    bool flag = false;
    bool flag2 = true;
    int numDegeneracies = 0;
    while((h!=comparisonHE || flag2) && (numDegeneracies < 2) ){
        // assert(numDegeneracies < 2);
        flag2=false;
        if(h->next()->twin()->vertex()->id() == v_new->id()){//degenerate face
            // if(h->twin()->vertex()->degree() <=2){return std::nullopt;}
            
            numDegeneracies++;
            FaceRef f1 = h->twin()->face();

            HalfedgeRef h_twin_prev = h->twin();
            while(h_twin_prev->next()!=h->twin()){
                h_twin_prev = h_twin_prev->next();
            }

            if(h->twin() == h->next()->twin()->next()){
                //delete both edges
                h->next()->twin()->face()->halfedge() = h->twin()->next();
                erase(h->next()->twin()); erase(h->next()); erase(h->next()->edge()); 
                //possibly if deg(other v) is 2 delete
            }
            else{
                currChanges[numDegeneracies-1].h = h;
                currChanges[numDegeneracies-1].h_twin_prev = h_twin_prev;

            }
            f1->halfedge() = h_twin_prev;
            
            toBeErased[numDegeneracies-1] = h;
            degenerateHE[numDegeneracies-1] = h;
            numDegenerateHE = numDegeneracies;
        }
        else{  // Non degenerate case
            if(h!=e->halfedge() && h!=e->halfedge()->twin()){
                h->face()->halfedge() = h;

                HalfedgeRef h_prev_prev = h;
                while(h_prev_prev->next()->next()!=h){
                    h_prev_prev = h_prev_prev->next();
                }

                if(h_prev_prev->next() == e->halfedge() || h_prev_prev->next() == e->halfedge()->twin()){
                    ChangesND lChangesND;
                    lChangesND.h = h;
                    lChangesND.h_prev_prev = h_prev_prev;
                    currChangesND.push_back(lChangesND);
                }
            }
            
        }
        h = h->twin()->next();

        if(!flag && h==comparisonHE){
            comparisonHE = e->halfedge()->twin()->next();
            h = comparisonHE;
            flag = true;
            flag2=true;
        }
    }

    for (int i=0; i<numDegeneracies; i++){
        HalfedgeRef h = toBeErased[i];
        erase(h->face()); erase(h->edge()); 
        erase(h->twin()); erase(h);
    }

    for (int i=0; i<numDegeneracies; i++){
        HalfedgeRef h = currChanges[i].h;
        HalfedgeRef h_twin_prev = currChanges[i].h_twin_prev;
        FaceRef f1 = h->twin()->face();

        h_twin_prev->next() = h->next();
        HalfedgeRef h_next_next = h->twin()->next();
        if(h_next_next == e->halfedge() || h_next_next == e->halfedge()->twin()){
            h_next_next = h_next_next->next();
        }
        h->next()->set_neighbors(h_next_next, h->next()->twin(), h->next()->vertex(), h->next()->edge(), f1);
        h->next()->vertex()->halfedge() = h->next();
        v_new->halfedge() = h->next()->twin();
    }

    for(auto lChangesND : currChangesND){
        lChangesND.h_prev_prev->next() = lChangesND.h;
        lChangesND.h->vertex()->halfedge() = lChangesND.h;
        v_new->halfedge() = lChangesND.h;
    }


    erase(e->halfedge()->twin()->vertex());
    erase(e->halfedge()->vertex());
    erase(e->halfedge()->twin());
    erase(e->halfedge());
    erase(e);
    
    return v_new;
}

/*
    This method should collapse the given face and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(Halfedge_Mesh::FaceRef f) {

    (void)f;
    return std::nullopt;
}

/*
    This method should flip the given edge and return an iterator to the
    flipped edge.
*/
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(Halfedge_Mesh::EdgeRef e) {
    HalfedgeRef h1 = e->halfedge(), h1_prev = h1;
    HalfedgeRef h2 = h1->twin(),h2_prev = h2;
    HalfedgeRef h1next = h1->next(), h2next = h2->next();
    while(h1_prev->next() != h1){
        h1_prev = h1_prev->next();
    }
    while(h2_prev->next() != h2){
        h2_prev = h2_prev->next();
    }
    //update vertices and face of h1,h2
    h1->vertex()->halfedge() = h2next;
    h2->vertex()->halfedge() = h1next;
    
    FaceRef f1 = h1->face();
    FaceRef f2 = h2->face();
    f1->halfedge() = h1;
    f2->halfedge() = h2;

    h1->set_neighbors(h1next->next(), h2, h2next->next()->vertex(), e, f1);
    h2->set_neighbors(h2next->next(), h1, h1next->next()->vertex(), e, f2);


    h1next->face() = f2;
    h2next->face() = f1;
    
    h1_prev->next() = h2next;
    h2_prev->next() = h1next;

    h1next->next() = h2;
    h2next->next() = h1;

    return e;
}

/*
    This method should split the given edge and return an iterator to the
    newly inserted vertex. The halfedge of this vertex should point along
    the edge that was split, rather than the new edges.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(Halfedge_Mesh::EdgeRef e) {
    HalfedgeRef h = e->halfedge();
    
    VertexRef v1 = h->vertex(), v2 = h->next()->vertex(),
              v3 = h->next()->next()->vertex();
    VertexRef v_new = new_vertex();
    v_new->pos = 0.5*v1->pos+ 0.5*v2->pos;
    HalfedgeRef h1 = new_halfedge(), h3 = new_halfedge(),
                h5 = new_halfedge(), h6 = new_halfedge(),
                h2 = new_halfedge(), h4 = new_halfedge();
    EdgeRef e1 = new_edge(), e2 = new_edge(), e3 = new_edge();
    FaceRef f1, f3;
    if(h->is_boundary()){//set boundary of faces
        f3 = new_face(h->next()->is_boundary());
        f1 = new_face(h->next()->next()->is_boundary());
    }

    f1->halfedge() = h->next()->next();
    f3->halfedge() = h->next();

    e1->halfedge() = h1;
    e2->halfedge() = h3;
    e3->halfedge() = h5;

    h1->set_neighbors(h5, h2, v1, e1, f1);
    h2->set_neighbors(h->twin()->next(), h1, v_new, e1, h->twin()->face()); //set face to be the virtual face
    h3->set_neighbors(h->next(),h4, v_new, e2, f3);
    h4->set_neighbors(h2, h3, v2, e2, h->twin()->face());

    HalfedgeRef h_twin_prev = h->twin()->next();
    while(h_twin_prev->next() != h->twin()){h_twin_prev = h_twin_prev->next();}
    h_twin_prev->set_neighbors(h4, h_twin_prev->twin(), h_twin_prev->vertex(), h_twin_prev->edge(), h_twin_prev->face());

    h5->set_neighbors(h->next()->next(), h6, v_new, e3, f1);
    h6->set_neighbors(h3, h5, v3, e3, f3);

    //edit old halfedges
    h->next()->next()->next() = h1;
    h->next()->next()->face() = f1;

    h->next()->next() = h6;
    h->next()->face() = f3;
    //edit v1 and v2
    v1->halfedge() = h1;
    v2->halfedge() = h->next();

    v_new->halfedge() = h3;

    if(!e->on_boundary()){
        VertexRef v4 = h->twin()->next()->next()->vertex();
        HalfedgeRef  h7 = new_halfedge(), h8 = new_halfedge();
        EdgeRef e4 = new_edge();
        FaceRef f2, f4;
        if(h->twin()->is_boundary()){//set boundary of faces
            f2 = new_face(h->twin()->next()->is_boundary());
            f4 = new_face(h->twin()->next()->next()->is_boundary());
        }
        f2->halfedge() = h->twin()->next();
        f4->halfedge() = h->twin()->next()->next();
        e4->halfedge() = h7;

        h1->set_neighbors(h5, h2, v1, e1, f1);
        h2->set_neighbors(h->twin()->next(), h1, v_new, e1, f2);
        h3->set_neighbors(h->next(), h4, v_new, e2, f3);
        h4->set_neighbors(h8, h3, v2, e2, f4);
        h7->set_neighbors(h2, h8, v4, e4, f2);
        h8->set_neighbors(h->twin()->next()->next(), h7, v_new, e4, f4);
        
        h->twin()->next()->next()->next() = h4;
        h->twin()->next()->next()->face() = f4;
        
        h->twin()->next()->next() = h7;
        h->twin()->next()->face() = f2;

        erase(h->twin()->face());
    }
    else{
        std::cout << "HE "<< h->id() << " is boundary" << std::endl;
    }
    erase(e); erase(h->face()); erase(h->twin()); erase(h);
    return v_new;
}

/* Note on the beveling process:

    Each of the bevel_vertex, bevel_edge, and bevel_face functions do not represent
    a full bevel operation. Instead, they should update the _connectivity_ of
    the mesh, _not_ the positions of newly created vertices. In fact, you should set
    the positions of new vertices to be exactly the same as wherever they "started from."

    When you click on a mesh element while in bevel mode, one of those three functions
    is called. But, because you may then adjust the distance/offset of the newly
    beveled face, we need another method of updating the positions of the new vertices.

    This is where bevel_vertex_positions, bevel_edge_positions, and
    bevel_face_positions come in: these functions are called repeatedly as you
    move your mouse, the position of which determins the normal and tangent offset
    parameters. These functions are also passed an array of the original vertex
    positions: for  bevel_vertex, it has one element, the original vertex position,
    for bevel_edge,  two for the two vertices, and for bevel_face, it has the original
    position of each vertex in halfedge order. You should use these positions, as well
    as the normal and tangent offset fields to assign positions to the new vertices.

    Finally, note that the normal and tangent offsets are not relative values - you
    should compute a particular new position from them, not a delta to apply.
*/

/*
    This method should replace the vertex v with a face, corresponding to
    a bevel operation. It should return the new face.  NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_vertex_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(Halfedge_Mesh::VertexRef v) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)v;
    return std::nullopt;
}

/*
    This method should replace the edge e with a face, corresponding to a
    bevel operation. It should return the new face. NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_edge_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(Halfedge_Mesh::EdgeRef e) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)e;
    return std::nullopt;
}

/*
    This method should replace the face f with an additional, inset face
    (and ring of faces around it), corresponding to a bevel operation. It
    should return the new face.  NOTE: This method is responsible for updating
    the *connectivity* of the mesh only---it does not need to update the vertex
    positions. These positions will be updated in
    Halfedge_Mesh::bevel_face_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_face(Halfedge_Mesh::FaceRef f) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."
    FaceRef f_new = new_face(f->is_boundary());
    
    HalfedgeRef h_in_prev = new_halfedge(), h_in_new, h_in_first = h_in_prev;
    f_new->halfedge() = h_in_prev;

    EdgeRef e_in_prev = new_edge();
    e_in_prev->halfedge() = h_in_prev;

    HalfedgeRef h = f->halfedge(), dummy;
    while(h->next()!=f->halfedge()){
        VertexRef v_prev = new_vertex();
        v_prev->pos = h->vertex()->pos;
        v_prev->halfedge() = h_in_prev;
        h_in_new = new_halfedge(); 
        h_in_prev->set_neighbors(h_in_new, dummy, v_prev, e_in_prev, f_new);
        h_in_prev = h_in_new;
        e_in_prev = new_edge();
        e_in_prev->halfedge() = h_in_new;
        h = h->next();
    }
    VertexRef v_prev = new_vertex();
    v_prev->pos = h->vertex()->pos;
    v_prev->halfedge() = h_in_prev;
    h_in_prev->set_neighbors(h_in_first, dummy, v_prev, e_in_prev, f_new);

    HalfedgeRef h_in = h_in_prev, h_in_orig = h_in;
    std::vector<HalfedgeRef> h_in_twins;
    int num_sides = 0;
    while(h_in->next()!=h_in_orig){
        FaceRef f_side = new_face(f->is_boundary());
        HalfedgeRef h_in_twin = new_halfedge();
        h_in->twin() = h_in_twin;
        h_in_twins.push_back(h_in_twin);
        h_in_twin->face() = f_side;
        f_side->halfedge() = h_in_twin;
        h_in = h_in->next();
        num_sides++;
    }
    FaceRef f_side = new_face(f->is_boundary());
    HalfedgeRef h_in_twin = new_halfedge();
    h_in->twin() = h_in_twin;
    h_in_twins.push_back(h_in_twin);
    h_in_twin->face() = f_side;
    f_side->halfedge() = h_in_twin;
    h_in = h_in->next();
    num_sides++;

    int i=0, idx;
    while(h_in->next()!=h_in_orig){
        idx = i-1<0 ? num_sides-1: i-1; 
        h_in->twin()->set_neighbors(h_in_twins[idx], h_in, h_in->next()->vertex(), h_in->edge(), h_in->twin()->face());
        i++;
        h_in = h_in->next();
    }
    idx = i-1<0 ? num_sides-1: i-1; 
    h_in->twin()->set_neighbors(h_in_twins[idx], h_in, h_in->next()->vertex(), h_in->edge(), h_in->twin()->face());
    h_in = h_in->next();

    //create edges connecting bevel face to original face
    EdgeRef e_prev = new_edge(), e_new, e_prev_orig = e_prev;
    HalfedgeRef h1_end = new_halfedge(), h1= h1_end,h2, hnext;
    e_prev->halfedge() = h1_end;

    while(h_in->next()!=h_in_orig){
        h2 = new_halfedge();
        h_in->twin()->next() = h2;
        h2->set_neighbors(h, h1, h_in->vertex(), e_prev, h_in->twin()->face());
        h1->twin() = h2; //h1 does not exist at the start
        h1 = new_halfedge();
        e_new = new_edge();
        e_new->halfedge() = h1;
        h1->set_neighbors(h_in->twin(), dummy, h->twin()->vertex(), e_new, h_in->twin()->face());
        
        e_prev = e_new;
        h->face() = h_in->twin()->face();
        h_in = h_in->next();
        hnext = h->next();
        h->next() = h1;
        h = hnext;
    }
    h2 = new_halfedge();
    h_in->twin()->next() = h2;
    h2->set_neighbors(h, h1, h_in->vertex(), e_prev, h_in->twin()->face());
    h1->twin() = h2;
    h1 = h1_end;
    e_prev_orig->halfedge() = h1;
    h1->set_neighbors(h_in->twin(), h_in_orig->twin()->next(), h->twin()->vertex(),e_prev_orig, h_in->twin()->face());

    h->face() = h_in->twin()->face();
    h->next() = h1;
    erase(f);
    return f_new;
}

/*
    Compute new vertex positions for the vertices of the beveled vertex.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the original vertex position and its associated outgoing edge
    to compute a new vertex position along the outgoing edge.
*/
void Halfedge_Mesh::bevel_vertex_positions(const std::vector<Vec3>& start_positions,
                                           Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
}

/*
    Compute new vertex positions for the vertices of the beveled edge.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the orig array) to compute an offset vertex position.

    Note that there is a 1-to-1 correspondence between halfedges in
    newHalfedges and vertex positions
    in orig.  So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vector3D pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_edge_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
}

/*
    Compute new vertex positions for the vertices of the beveled face.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the start_positions array) to compute an offset vertex
    position.

    Note that there is a 1-to-1 correspondence between halfedges in
    new_halfedges and vertex positions
    in orig. So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vec3 pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_face_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset,
                                         float normal_offset) {

    if(flip_orientation) normal_offset = -normal_offset;
    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());
    
    Vec3 avg_pos = Vec3();
    for(size_t i = 0; i < start_positions.size(); i++){
            avg_pos += start_positions[i]; 
    }
    avg_pos /= start_positions.size();

    float alpha = std::max(1e-2, 2.0*tangent_offset/3 + 1.0);
    //std::cout<<"tangent_offset="<<tangent_offset<<"alpha ="<<alpha<<std::endl;
    Vec3 shift = face->normal()*normal_offset;

    for(size_t i = 0; i < new_halfedges.size(); i++){
        new_halfedges[i]->vertex()->pos = start_positions[i]*alpha+avg_pos*(1-alpha) + shift;
    }
    if(alpha<=1e-2){
        collapse_face(face);
    }
    


    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
    (void)normal_offset;
}

/*
    Splits all non-triangular faces into triangles.
*/
void Halfedge_Mesh::triangulate() {
    for(FaceRef f = faces_begin(); f != faces_end(); f++) {
        if(f->is_boundary()){continue;}

        HalfedgeRef h = f->halfedge();
        HalfedgeRef  h0 = h, h1 = h->next();
        while(h1->next()->next() != h){
            HalfedgeRef hnew = new_halfedge(), hnew_twin = new_halfedge(), h1next = h1->next(); // hnew is the halfedge inside the new triangle
            EdgeRef enew = new_edge();
            FaceRef fnew = new_face(f->is_boundary());
            
            hnew->set_neighbors(h0, hnew_twin, h1->twin()->vertex(), enew, fnew);
            hnew_twin->set_neighbors(h1->next(), hnew, h0->vertex(), enew, f);

            h0->set_neighbors(h1, h0->twin(), h0->vertex(), h0->edge(), fnew);
            h1->set_neighbors(hnew, h1->twin(), h1->vertex(), h1->edge(), fnew);

            enew->halfedge() = hnew;
            fnew->halfedge() = hnew;

            h0 = hnew_twin;
            h1 = h1next;
        }
        f->halfedge() = h1;
        h1->next()->next() = h0;
    }    
}

/* Note on the quad subdivision process:

        Unlike the local mesh operations (like bevel or edge flip), we will perform
        subdivision by splitting *all* faces into quads "simultaneously."  Rather
        than operating directly on the halfedge data structure (which as you've
        seen is quite difficult to maintain!) we are going to do something a bit nicer:
           1. Create a raw list of vertex positions and faces (rather than a full-
              blown halfedge mesh).
           2. Build a new halfedge mesh from these lists, replacing the old one.
        Sometimes rebuilding a data structure from scratch is simpler (and even
        more efficient) than incrementally modifying the existing one.  These steps are
        detailed below.

  Step I: Compute the vertex positions for the subdivided mesh.
        Here we're going to do something a little bit strange: since we will
        have one vertex in the subdivided mesh for each vertex, edge, and face in
        the original mesh, we can nicely store the new vertex *positions* as
        attributes on vertices, edges, and faces of the original mesh. These positions
        can then be conveniently copied into the new, subdivided mesh.
        This is what you will implement in linear_subdivide_positions() and
        catmullclark_subdivide_positions().

  Steps II-IV are provided (see Halfedge_Mesh::subdivide()), but are still detailed
  here:

  Step II: Assign a unique index (starting at 0) to each vertex, edge, and
        face in the original mesh. These indices will be the indices of the
        vertices in the new (subdivided mesh).  They do not have to be assigned
        in any particular order, so long as no index is shared by more than one
        mesh element, and the total number of indices is equal to V+E+F, i.e.,
        the total number of vertices plus edges plus faces in the original mesh.
        Basically we just need a one-to-one mapping between original mesh elements
        and subdivided mesh vertices.

  Step III: Build a list of quads in the new (subdivided) mesh, as tuples of
        the element indices defined above. In other words, each new quad should be
        of the form (i,j,k,l), where i,j,k and l are four of the indices stored on
        our original mesh elements.  Note that it is essential to get the orientation
        right here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces
        should circulate in the same direction as old faces (think about the right-hand
        rule).

  Step IV: Pass the list of vertices and quads to a routine that clears
        the internal data for this halfedge mesh, and builds new halfedge data from
        scratch, using the two lists.
*/

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    simple linear interpolation, e.g., the edge midpoints and face
    centroids.
*/
void Halfedge_Mesh::linear_subdivide_positions() {

    // For each vertex, assign Vertex::new_pos to
    // its original position, Vertex::pos.
    for(VertexRef v = this->vertices_begin(); v != this->vertices_end(); v++) {
        v->setNewPos(v->pos);
    }

    // For each edge, assign the midpoint of the two original
    // positions to Edge::new_pos.
    for(EdgeRef e = this->edges_begin(); e != this->edges_end(); e++) {
        e->setNewPos(0.5*e->halfedge()->vertex()->pos + 0.5*e->halfedge()->twin()->vertex()->pos);
    }

    // For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::new_pos. Note
    // that in general, NOT all faces will be triangles!
    for(FaceRef f = this->faces_begin(); f != this->faces_end(); f++) {
        std::vector<Vec3> vertices;
        HalfedgeRef h = f->halfedge();
        while(h->next() != f->halfedge()){
            vertices.push_back(h->vertex()->pos);
            h=h->next();
        }
        vertices.push_back(h->vertex()->pos);
        
        Vec3 newPos(0,0,0);
        for(auto vertex: vertices){
            newPos += vertex / vertices.size();
        }
        f->setNewPos(newPos);
    }
}

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    the Catmull-Clark rules for subdivision.

    Note: this will only be called on meshes without boundary
*/
void Halfedge_Mesh::catmullclark_subdivide_positions() {

    // The implementation for this routine should be
    // a lot like Halfedge_Mesh:linear_subdivide_positions:(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules. (These rules are outlined in the Developer Manual.)

    // Faces
    for(FaceRef f = this->faces_begin(); f != this->faces_end(); f++) {
        std::vector<Vec3> vertices;
        HalfedgeRef h = f->halfedge();
        while(h->next() != f->halfedge()){
            vertices.push_back(h->vertex()->pos);
            h=h->next();
        }
        vertices.push_back(h->vertex()->pos);
        
        Vec3 newPos(0,0,0);
        for(auto vertex: vertices){
            newPos += vertex / vertices.size();
        }
        f->setNewPos(newPos);
    }

    // Edges
    for(EdgeRef e = this->edges_begin(); e != this->edges_end(); e++) {
        std::vector<Vec3> vertices;
        vertices.push_back(e->halfedge()->vertex()->pos);
        vertices.push_back(e->halfedge()->twin()->vertex()->pos);
        vertices.push_back(e->halfedge()->face()->getNewPos());
        vertices.push_back(e->halfedge()->twin()->face()->getNewPos());
        Vec3 newPos(0,0,0);
        for(auto vertex: vertices){
            newPos += vertex / vertices.size();
        }
        e->setNewPos(newPos);
    }

    // Vertices
    for(VertexRef v = this->vertices_begin(); v != this->vertices_end(); v++) {
        std::vector<Vec3> edgeMidPoints, faceNewPos;
        int degree=0;

        HalfedgeRef h = v->halfedge();
        while(h->twin()->next() != v->halfedge()){
            edgeMidPoints.push_back(0.5*(h->vertex()->pos + h->twin()->vertex()->pos));
            faceNewPos.push_back(h->face()->getNewPos());
            degree++;
            h = h->twin()->next(); 
        }
        edgeMidPoints.push_back(0.5*(h->vertex()->pos + h->twin()->vertex()->pos));
        faceNewPos.push_back(h->face()->getNewPos());
        degree++;

        Vec3 edgeNewPosAvg(0,0,0);
        for(auto newPos: edgeMidPoints){
            edgeNewPosAvg += newPos / edgeMidPoints.size();
        }

        Vec3 faceNewPosAvg(0,0,0);
        for(auto newPos: faceNewPos){
            faceNewPosAvg += newPos / faceNewPos.size();
        }

        Vec3 newPos = (faceNewPosAvg + 2*edgeNewPosAvg + (degree-3)*v->pos)/degree;
        v->setNewPos(newPos);
    }
}

/*
        This routine should increase the number of triangles in the mesh
        using Loop subdivision. Note: this is will only be called on triangle meshes.
*/
void Halfedge_Mesh::loop_subdivide() {

    // Compute new positions for all the vertices in the input mesh, using
    // the Loop subdivision rule, and store them in Vertex::new_pos.
    // -> At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh. Use Vertex::is_new for this.
    // -> Next, compute the updated vertex positions associated with edges, and
    //    store it in Edge::new_pos.
    // -> Next, we're going to split every edge in the mesh, in any order.  For
    //    future reference, we're also going to store some information about which
    //    subdivided edges come from splitting an edge in the original mesh, and
    //    which edges are new, by setting the flat Edge::is_new. Note that in this
    //    loop, we only want to iterate over edges of the original mesh.
    //    Otherwise, we'll end up splitting edges that we just split (and the
    //    loop will never end!)
    // -> Now flip any new edge that connects an old and new vertex.
    // -> Finally, copy the new vertex positions into final Vertex::pos.

    // Each vertex and edge of the original surface can be associated with a
    // vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to
    // *first* compute the new positions
    // using the connectivity of the original (coarse) mesh; navigating this mesh
    // will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.  We
    // will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.

    // Compute updated positions for all the vertices in the original mesh, using
    // the Loop subdivision rule.

    // Next, compute the updated vertex positions associated with edges.

    // Next, we're going to split every edge in the mesh, in any order. For
    // future reference, we're also going to store some information about which
    // subdivided edges come from splitting an edge in the original mesh, and
    // which edges are new.
    // In this loop, we only want to iterate over edges of the original
    // mesh---otherwise, we'll end up splitting edges that we just split (and
    // the loop will never end!)

    // Finally, flip any new edge that connects an old and new vertex.

    // Copy the updated vertex positions to the subdivided mesh.
}

/*
    Isotropic remeshing. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if this is not a triangle mesh)
*/
bool Halfedge_Mesh::isotropic_remesh() {

    // Compute the mean edge length.
    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    return false;
}

/* Helper type for quadric simplification */
struct Edge_Record {
    Edge_Record() {
    }
    Edge_Record(std::unordered_map<Halfedge_Mesh::VertexRef, Mat4>& vertex_quadrics,
                Halfedge_Mesh::EdgeRef e)
        : edge(e) {

        // Compute the combined quadric from the edge endpoints.
        // -> Build the 3x3 linear system whose solution minimizes the quadric error
        //    associated with these two endpoints.
        // -> Use this system to solve for the optimal position, and store it in
        //    Edge_Record::optimal.
        // -> Also store the cost associated with collapsing this edge in
        //    Edge_Record::cost.

        Mat4 edgeK = vertex_quadrics[e->halfedge()->vertex()] + vertex_quadrics[e->halfedge()->twin()->vertex()];
        Vec3 b(-edgeK[0][3], -edgeK[1][3], -edgeK[2][3]);
        edgeK[0][3]=0;edgeK[3][0]=0;edgeK[1][3]=0;edgeK[3][1]=0;edgeK[2][3]=0;edgeK[3][2]=0;edgeK[3][3]=1;
        this->optimal = edgeK.inverse() * b;
        
        // Vec3 temp = this->optimal * (edgeK * this->optimal);
        // this->cost = temp.x+temp.y+temp.z;
        this->cost = dot(this->optimal, (edgeK * this->optimal));
        // std::cout << "EdgeRecord edge " << edge->id() << " cost " << cost << std::endl;

    }
    Halfedge_Mesh::EdgeRef edge;
    Vec3 optimal;
    float cost;
};

/* Comparison operator for Edge_Records so std::set will properly order them */
bool operator<(const Edge_Record& r1, const Edge_Record& r2) {
    if(r1.cost != r2.cost) {
        return r1.cost < r2.cost;
    }
    Halfedge_Mesh::EdgeRef e1 = r1.edge;
    Halfedge_Mesh::EdgeRef e2 = r2.edge;
    return &*e1 < &*e2;
}

/** Helper type for quadric simplification
 *
 * A PQueue is a minimum-priority queue that
 * allows elements to be both inserted and removed from the
 * queue.  Together, one can easily change the priority of
 * an item by removing it, and re-inserting the same item
 * but with a different priority.  A priority queue, for
 * those who don't remember or haven't seen it before, is a
 * data structure that always keeps track of the item with
 * the smallest priority or "score," even as new elements
 * are inserted and removed.  Priority queues are often an
 * essential component of greedy algorithms, where one wants
 * to iteratively operate on the current "best" element.
 *
 * PQueue is templated on the type T of the object
 * being queued.  For this reason, T must define a comparison
 * operator of the form
 *
 *    bool operator<( const T& t1, const T& t2 )
 *
 * which returns true if and only if t1 is considered to have a
 * lower priority than t2.
 *
 * Basic use of a PQueue might look
 * something like this:
 *
 *    // initialize an empty queue
 *    PQueue<myItemType> queue;
 *
 *    // add some items (which we assume have been created
 *    // elsewhere, each of which has its priority stored as
 *    // some kind of internal member variable)
 *    queue.insert( item1 );
 *    queue.insert( item2 );
 *    queue.insert( item3 );
 *
 *    // get the highest priority item currently in the queue
 *    myItemType highestPriorityItem = queue.top();
 *
 *    // remove the highest priority item, automatically
 *    // promoting the next-highest priority item to the top
 *    queue.pop();
 *
 *    myItemType nextHighestPriorityItem = queue.top();
 *
 *    // Etc.
 *
 *    // We can also remove an item, making sure it is no
 *    // longer in the queue (note that this item may already
 *    // have been removed, if it was the 1st or 2nd-highest
 *    // priority item!)
 *    queue.remove( item2 );
 *
 */
template<class T> struct PQueue {
    void insert(const T& item) {
        queue.insert(item);
    }
    void remove(const T& item) {
        if(queue.find(item) != queue.end()) {
            queue.erase(item);
        }
    }
    const T& top(void) const {
        return *(queue.begin());
    }
    void pop(void) {
        queue.erase(queue.begin());
    }
    size_t size() {
        return queue.size();
    }

    std::set<T> queue;
};

/*
    Mesh simplification. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if you can't simplify the mesh any
    further without destroying it.)
*/
bool Halfedge_Mesh::simplify() {

    std::unordered_map<VertexRef, Mat4> vertex_quadrics;
    std::unordered_map<FaceRef, Mat4> face_quadrics;
    std::unordered_map<EdgeRef, Edge_Record> edge_records;
    PQueue<Edge_Record> edge_queue;

    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in face_quadrics
    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in vertex_quadrics
    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an Edge_Record for each edge and sticking it in the
    //    queue. You may want to use the above PQueue<Edge_Record> for this.
    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    /* Compute Face Quadric Matrix K */
    std::cout << "/* Compute Face Quadric Matrix K */" << std::endl;
    for(FaceRef f = faces_begin(); f != faces_end(); f++) {        
        Vec3 normal = f->normal();
        Vec3 vertexPos = f->halfedge()->vertex()->pos;
        float d = -1*(normal.x*vertexPos.x + normal.y*vertexPos.y + normal.z*vertexPos.z);
        Vec4 v(normal.x, normal.y, normal.z, d);
        Mat4 K = outer(v, v);
        face_quadrics[f] = K;
    }

    /* Compute Vertex Quadric for each vertex */
    std::cout << "/* Compute Vertex Quadric for each vertex */" << std::endl;
    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {        
        /* Iterate over all faces touching the vertex */
        HalfedgeRef h = v->halfedge();
        Mat4 K = Mat4::Zero;
        while(h->twin()->next() != v->halfedge()){
            if(!h->face()->is_boundary()){
                K += face_quadrics[h->face()];
            }
            h=h->twin()->next();
        }
        if(!h->face()->is_boundary()){
            K += face_quadrics[h->face()];
        }
        vertex_quadrics[v] = K;
    }

    /* EdgeRecord for each edge */
    std::cout << "/* EdgeRecord for each edge */" << std::endl;
    for(EdgeRef e = edges_begin(); e != edges_end(); e++) {        
        Edge_Record newER(vertex_quadrics, e);
        edge_records[e] = newER;
        edge_queue.insert(newER);
    }

    /* Collapse best edge repeatedly */
    std::cout << "/* Collapse best edge repeatedly */" << std::endl;
    const unsigned int targetNumEdges = this->edges.size()/4;
    int numIterations = 0;
    while(this->edges.size() > targetNumEdges){
        std::cout << "Inside While: NumEdges="<<this->edges.size() << " Target Num Edges: " << targetNumEdges << " NumIterations: "<< numIterations << std::endl;
        // Get the cheapest edge from the queue.
        Edge_Record bestER = edge_queue.top();
        edge_queue.pop();

        EdgeRef bestEdge = bestER.edge;
        VertexRef v0 = bestEdge->halfedge()->vertex(), v1 = bestEdge->halfedge()->twin()->vertex();
        Mat4 newQuadric = vertex_quadrics[v0] + vertex_quadrics[v1];

        // Remove any edge touching either of its endpoints from the queue.
        std::cout << "Remove any edge touching either of its endpoints from the queue.: V0" << std::endl;
        HalfedgeRef h = v0->halfedge();
        while(h->twin()->next() != v0->halfedge()){
            edge_queue.remove(edge_records[h->edge()]);
            h = h->twin()->next();
        }
        edge_queue.remove(edge_records[h->edge()]);

        std::cout << "Remove any edge touching either of its endpoints from the queue.: V1" << std::endl;
        h = v1->halfedge();
        while(h->twin()->next() != v1->halfedge()){
            edge_queue.remove(edge_records[h->edge()]);
            h = h->twin()->next();
        }
        edge_queue.remove(edge_records[h->edge()]);

        // Collapse the edge.
        // assert(edge_queue.queue.find(bestER) == edge_queue.queue.end());
        std::cout << "Collapsing edge "<<bestEdge->id() <<std::endl;
        std::optional<VertexRef> newVertexOpt = collapse_edge_erase(bestEdge);
        if(newVertexOpt == std::nullopt){
            std::cout << "newVertexOpt == std::nullopt, Cannot Collapse Edge " << bestEdge->id() << std::endl;
            continue;
        } 
        for(int i=0; i<numDegenerateHE; i++){
            edge_queue.remove(edge_records[degenerateHE[i]->edge()]);
        }
        VertexRef newVertex = newVertexOpt.value();

        // Set the quadric of the new vertex to the quadric computed in Step 3.
        std::cout << "New vertex created after edge collapse " << newVertex->id() << std::endl;
        vertex_quadrics[newVertex] = newQuadric;

        // Insert any edge touching the new vertex into the queue, creating new edge records for each of them.
        std::cout << "Insert any edge touching the new vertex into the queue, creating new edge records for each of them." << std::endl;
        h = newVertex->halfedge();
        // int hId = h->id(); (void) hId;
        while(h->twin()->next() != newVertex->halfedge()){            
            Edge_Record newER(vertex_quadrics, h->edge());
            edge_records[h->edge()] = newER;
            edge_queue.insert(newER);
            h = h->twin()->next();
            // hId = h->id(); (void) hId;
        }
        // edge_queue.remove(edge_records[h->edge()]);
        Edge_Record newER(vertex_quadrics, h->edge());
        edge_records[h->edge()] = newER;
        edge_queue.insert(newER);
        
        std::cout << "Exiting while Loop" << std::endl;
        numIterations ++;
        // if(numIterations ==100) break;
        // break;

        auto valResult = validate();
        if(valResult!=std::nullopt){
            std::pair<Halfedge_Mesh::ElementRef, std::string> object = valResult.value();
            HalfedgeRef elem = std::get<HalfedgeRef>(object.first);
            std::cout << elem->id() << " " << object.second <<std::endl;
            break;
        }
    }

    /* Collapse best edge repeatedly */
    std::cout << "Exit Simplify"<<std::endl;

    return true;
}
