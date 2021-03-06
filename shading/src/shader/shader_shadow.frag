//
// Parameters that control fragment shader behavior. Different materials
// will set these flags to true/false for different looks
//

uniform bool useTextureMapping;     // true if basic texture mapping (diffuse) should be used
uniform bool useNormalMapping;      // true if normal mapping should be used
uniform bool useEnvironmentMapping; // true if environment mapping should be used
uniform bool useMirrorBRDF;         // true if mirror brdf should be used (default: phong)

//
// texture maps
//

uniform sampler2D diffuseTextureSampler;

// TODO CS248 Part 3: Normal Mapping
uniform sampler2D normalTextureSampler;
// TODO CS248 Part 4: Environment Mapping
uniform sampler2D environmentSampler;
uniform sampler2DArray shadowSampler;

//
// lighting environment definition. Scenes may contain directional
// and point light sources, as well as an environment map
//

#define MAX_NUM_LIGHTS 10
uniform int  num_directional_lights;
uniform vec3 directional_light_vectors[MAX_NUM_LIGHTS];

uniform int  num_point_lights;
uniform vec3 point_light_positions[MAX_NUM_LIGHTS];

uniform int   num_spot_lights;
uniform vec3  spot_light_positions[MAX_NUM_LIGHTS];
uniform vec3  spot_light_directions[MAX_NUM_LIGHTS];
uniform vec3  spot_light_intensities[MAX_NUM_LIGHTS];
uniform float spot_light_angles[MAX_NUM_LIGHTS];


//
// material-specific uniforms
//

// parameters to Phong BRDF
uniform float spec_exp;

// values that are varying per fragment (computed by the vertex shader)

in vec3 position;     // surface position
in vec2 texcoord;     // surface texcoord (uv)
in vec3 normal;       // surface normal
in vec3 dir2camera;   // vector from surface point to camera
in mat3 tan2world;    // tangent space to world space transform
in vec3 vertex_diffuse_color; // surface color
in vec4 lightspace_position[MAX_NUM_LIGHTS];

out vec4 fragColor;

#define PI 3.14159265358979323846


//
// Simple diffuse brdf
//
// L -- direction to light
// N -- surface normal at point being shaded
//
vec3 Diffuse_BRDF(vec3 L, vec3 N, vec3 diffuseColor) {
    return diffuseColor * max(dot(N, L), 0.);
}

//
// Phong_BRDF --
//
// Evaluate phong reflectance model according to the given parameters
// L -- direction to light
// V -- direction to camera (view direction)
// N -- surface normal at point being shaded
//
vec3 Phong_BRDF(vec3 L, vec3 V, vec3 N, vec3 diffuse_color, vec3 specular_color, float specular_exponent)
{
    // TODO CS248 Part 2: Phong Reflectance
    // Implement diffuse and specular terms of the Phong
    // reflectance model here.
    vec3 L_norm = normalize(L);
    vec3 N_norm = normalize(N);
    vec3 V_norm = normalize(V);
    vec3 R_norm = normalize(2*dot(L_norm,N_norm)*N_norm-L_norm);
    float diffused_comp = max(dot(L_norm,N_norm),0.0);
    float spec_comp = pow(max(dot(R_norm,V_norm),0.0),specular_exponent);

    return diffused_comp*diffuse_color+spec_comp*specular_color;

}

//
// Blinn-Phong_BRDF --
//
// Evaluate phong reflectance model according to the given parameters
// L -- direction to light
// H -- direction of Halfway Vector
// N -- surface normal at point being shaded
//
vec3 Blinn_Phong_BRDF(vec3 L, vec3 H, vec3 N, vec3 diffuse_color, vec3 specular_color, float specular_exponent)
{
    vec3 L_norm = normalize(L);
    vec3 H_norm = normalize(L);
    vec3 N_norm = normalize(N);
    float diffused_comp = max(dot(L_norm,N_norm),0.0);
    float spec_comp = pow(max(dot(N_norm,H_norm),0.0),specular_exponent);

    return diffused_comp*diffuse_color+spec_comp*specular_color;
}

//
// Oren Nayar BRDF --
// http://www.cs.cmu.edu/afs/cs/academic/class/16823-s16/www/pdfs/appearance-modeling-5.pdf
// Evaluate reflectance model according to the given parameters
// L -- direction to light
// V -- direction to camera (view direction)
// N -- surface normal at point being shaded
// roughness -- the std. deviation of the normals directions of the surfece facets
vec3 Oren_Nayar_BRDF(vec3 L, vec3 V, vec3 N, float roughness, vec3 diffuse_color){
    vec3 L_norm = normalize(L);
    vec3 N_norm = normalize(N);
    vec3 V_norm = normalize(V);
    float cosI = dot(N_norm,L_norm), cosR = dot(N_norm,V_norm);

    float sigma2 = roughness*roughness;
    float A = 1 - 0.5*sigma2/(sigma2+0.33), B = 0.45*sigma2/(sigma2+0.09);
    
    float cosPhiSinAtanB = (dot(L_norm,V_norm) - cosI*cosR)/max(cosI,cosR);

    return max(dot(L_norm,N_norm),0.0)*diffuse_color*(A+B*max(0.0,cosPhiSinAtanB));
}

//
// Lommel Seeligger BRDF --
// https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Book%3A_Planetary_Photometry_(Tatum_and_Fairbairn)/03%3A_A_Brief_History_of_the_Lommel-Seeliger_Law/3.01%3A_A_Brief_History_of_the_Lommel-Seeliger_Law
//
// Evaluate  reflectance model according to the given parameters
// L -- direction to light
// V -- direction to camera (view direction)
// N -- surface normal at point being shaded
//
vec3 LS_BRDF(vec3 L, vec3 V, vec3 N, vec3 diffuse_color)
{
    vec3 L_norm = normalize(L);
    vec3 N_norm = normalize(N);
    vec3 V_norm = normalize(V);
    
    float u0 = max(dot(L_norm,N_norm),0.0);
    float u = max(dot(V_norm,N_norm),0.0);
    float brdf = 1.2* u0 / ((u + u0)); // Proportionality constant changed for more pleasant image

    return brdf*diffuse_color;

}

//
// Ward BRDF --
// https://www.graphics.cornell.edu/~bjw/wardnotes.pdf
//
// Evaluate  reflectance model according to the given parameters
// L -- direction to light
// V -- direction to camera (view direction)
// N -- surface normal at point being shaded
//
vec3 Ward_BRDF(vec3 L, vec3 V, vec3 N, vec3 diffuse_color, vec3 specular_color)
{
    vec3 L_norm = normalize(L);
    vec3 N_norm = normalize(N);
    vec3 V_norm = normalize(V);
    float diffused_comp = max(dot(L_norm,N_norm),0.0);
    

    const float rho_s = 1, alpha = 1, k = 0.3; // Assuming BRDF is isotropic: alpha_x = alpha_y
    vec3 halfDir = normalize(L+V);
    float cosThetaL = dot(L_norm, N_norm);
    float cosThetaV = dot(V_norm, N_norm);
    float cosThetaH = dot(halfDir, N_norm);
    float tanSqThetaH = -1 + 1.0/(cosThetaH * cosThetaH);
    float spec_comp = k* rho_s * pow(2.718, -1 * tanSqThetaH) / (alpha * alpha * sqrt(cosThetaL * cosThetaV));

    return diffused_comp*diffuse_color+spec_comp*specular_color;

}

//
// SampleEnvironmentMap -- returns incoming radiance from specified direction
//
// D -- world space direction (outward from scene) from which to sample radiance
// 
vec3 SampleEnvironmentMap(vec3 D)
{    
    // TODO CS248 Part 4: Environment Mapping
    // sample environment map in direction D.  This requires
    // converting D into spherical coordinates where Y is the polar direction
    // (warning: in our scene, theta is angle with Y axis, which differs from
    // typical convention in physics)
    //
    // Tips:
    //
    // (1) See GLSL documentation of acos(x) and atan(x, y)
    //
    // (2) atan() returns an angle in the range -PI to PI, so you'll have to
    //     convert negative values to the range 0 - 2PI
    //
    // (3) How do you convert theta and phi to normalized texture
    //     coordinates in the domain [0,1]^2?

    // vec3 normalized = normalize(D);
    float theta =  acos(D.y/(length(D)));

    float temp = atan(D.x, D.z)/(PI) + 2; // Range [1, 3)
    float phi = 2*PI*(temp/2 - int(temp/2)); // Range [0, 2*PI)

    // Between [0,1]^2
    phi/=2*PI;
    theta/=PI;
    
    vec3 color = texture(environmentSampler, vec2(phi, theta)).rgb;

    return color;    
}

//
// Fragment shader main entry point
//
void main(void)
{
    //////////////////////////////////////////////////////////////////////////
	// Phase 1: Pattern generation. Compute parameters to BRDF 
    //////////////////////////////////////////////////////////////////////////
    
	vec3 diffuseColor = vec3(1.0, 1.0, 1.0);
    vec3 specularColor = vec3(1.0, 1.0, 1.0);
    float specularExponent = spec_exp;

    if (useTextureMapping) {
        diffuseColor = texture(diffuseTextureSampler, texcoord).rgb;
    } else {
        diffuseColor = vertex_diffuse_color;
    }



    // perform normal map lookup if required
    vec3 N = vec3(0);
    if (useNormalMapping) {
       // TODO: CS248 Part 3: Normal Mapping:
       // use tan2World in the normal map to compute the
       // world space normal baaed on the normal map.

       // Note that values from the texture should be scaled by 2 and biased
       // by negative -1 to covert positive values from the texture fetch, which
       // lie in the range (0-1), to the range (-1,1).
       //
       // In other words:   tangent_space_normal = texture_value * 2.0 - 1.0;

       // replace this line with your implementation
       N = texture(normalTextureSampler, texcoord).rgb;
       N = 2*N-1.0;
       N = tan2world*N;

    } else {
       N = normalize(normal);
    }

    vec3 V = normalize(dir2camera);
    vec3 Lo = vec3(0.1 * diffuseColor);   // this is ambient

    /////////////////////////////////////////////////////////////////////////
    // Phase 2: Evaluate lighting and surface BRDF 
    /////////////////////////////////////////////////////////////////////////

    if (useMirrorBRDF) {
        //
        // TODO: CS248 Part 4: Environment Mapping:
        // compute perfect mirror reflection direction here.
        // You'll also need to implement environment map sampling in SampleEnvironmentMap()
        //
        vec3 normalizedN = normalize(normal);
        vec3 alongNormal = dot(dir2camera, normalizedN) * normalizedN;
        vec3 awayFromNormal = dir2camera - alongNormal;
        vec3 R = alongNormal - awayFromNormal;
        //

        // sample environment map
        vec3 envColor = SampleEnvironmentMap(R);
        
        // this is a perfect mirror material, so we'll just return the light incident
        // from the reflection direction
        fragColor = vec4(envColor, 1);
        return;
    }

	// for simplicity, assume all lights (other than spot lights) have unit magnitude
	float light_magnitude = 1.0;

	// for all directional lights
	for (int i = 0; i < num_directional_lights; ++i) {
	    vec3 L = normalize(-directional_light_vectors[i]);
		vec3 brdf_color = Phong_BRDF(L, V, N, diffuseColor, specularColor, specularExponent);
        //Lommel Seeligger Reflectance
        // vec3 brdf_color = LS_BRDF(L, V, N, diffuseColor);

        //Ward BRDF
        // vec3 brdf_color = Ward_BRDF(L, V, N, diffuseColor, specular_color);
        
        //Blinn-Phong Reflectance
        // vec3 H = normalize(L+V);
        // vec3 brdf_color = Blinn_Phong_BRDF(L, H, N, diffuseColor, specularColor, specularExponent);
	    
        //Oren-Nayar Reflectance
        // float roughness=0.5;
        // vec3 brdf_color = Oren_Nayar_BRDF(L, V, N, roughness, diffuseColor);
	    Lo += light_magnitude * brdf_color;
    }

    // for all point lights
    for (int i = 0; i < num_point_lights; ++i) {
		vec3 light_vector = point_light_positions[i] - position;
        vec3 L = normalize(light_vector);
        float distance = length(light_vector);
        vec3 brdf_color = Phong_BRDF(L, V, N, diffuseColor, specularColor, specularExponent);
        //Lommel Seeligger Reflectance
        // vec3 brdf_color = LS_BRDF(L, V, N, diffuseColor);

        //Ward BRDF
        // vec3 brdf_color = Ward_BRDF(L, V, N, diffuseColor, specular_color);
        
        //Blinn-Phong Reflectance
        // vec3 H = normalize(L+V);
        // vec3 brdf_color = Blinn_Phong_BRDF(L, H, N, diffuseColor, specularColor, specularExponent);
	    
        //Oren-Nayar Reflectance
        // float roughness=0.5;
        // vec3 brdf_color = Oren_Nayar_BRDF(L, V, N, roughness, diffuseColor);
        float falloff = 1.0 / (0.01 + distance * distance);
        Lo += light_magnitude * falloff * brdf_color;
    }

    // for all spot lights
	for (int i = 0; i < num_spot_lights; ++i) {
    
        vec3 intensity = spot_light_intensities[i];   // intensity of light: this is intensity in RGB
        vec3 light_pos = spot_light_positions[i];     // location of spotlight
        float cone_angle = spot_light_angles[i];      // spotlight falls off to zero in directions whose
                                                      // angle from the light direction is grester than
                                                      // cone angle. Caution: this value is in units of degrees!

        vec3 dir_to_surface = position - light_pos;
        float angle = acos(dot(normalize(dir_to_surface), spot_light_directions[i])) * 180.0 / PI;

        // TODO CS248 Part 5.1: Spotlight Attenuation: compute the attenuation of the spotlight due to two factors:
        // (1) distance from the spot light (D^2 falloff)
        // (2) attentuation due to being outside the spotlight's cone 
        //
        // Here is a description of what to compute:
        //
        // 1. Modulate intensity by a factor of 1/D^2, where D is the distance from the
        //    spotlight to the current surface point.  For robustness, it's common to use 1/(1 + D^2)
        //    to never multiply by a value greather than 1.
        //
        // 2. Modulate the resulting intensity based on whether the surface point is in the cone of
        //    illumination.  To achieve a smooth falloff, consider the following rules
        //    
        //    -- Intensity should be zero if angle between the spotlight direction and the vector from
        //       the light position to the surface point is greater than (1.0 + SMOOTHING) * cone_angle
        //
        //    -- Intensity should not be further attentuated if the angle is less than (1.0 - SMOOTHING) * cone_angle
        //
        //    -- For all other angles between these extremes, interpolate linearly from unattenuated
        //       to zero intensity. 
        //
        //    -- The reference solution uses SMOOTHING = 0.1, so 20% of the spotlight region is the smoothly
        //       facing out area.  Smaller values of SMOOTHING will create hard spotlights.

        // CS248: remove this once you perform proper attenuation computations
        float D = length(dir_to_surface);
        intensity*=1/(1+D*D);

        float SMOOTHING = 0.1;
        if(angle > (1.0 + SMOOTHING)* cone_angle){
            intensity = vec3(0,0,0);
        }
        else if(angle > (1.0-SMOOTHING)*cone_angle){
            float t = ((1.0 + SMOOTHING)* cone_angle - angle) / (2*SMOOTHING*cone_angle);
            intensity *= t;
        }


        // Render Shadows for all spot lights
        // TODO CS248 Part 5.2: Shadow Mapping: comute shadowing for spotlight i here 
        vec4 position_shadowlight = lightspace_position[i];
        vec3 lightspace3D = position_shadowlight.xyz / position_shadowlight.w;
        vec2 shadow_uv = lightspace3D.xy;
        float light2PointDistance = lightspace3D.z;

        float pcf_step_size = 256;
        int num_ilum = 0;
        for (int j=-2; j<=2; j++) {
            for (int k=-2; k<=2; k++) {
                vec2 offset = vec2(j,k) / pcf_step_size;
                float closestDistance = texture(shadowSampler, vec3(shadow_uv+offset, i)).x;
                if(closestDistance >= light2PointDistance- 0.005){
                    num_ilum++;
                }
            }
        }
        float frac = float(num_ilum)/25.0;
        intensity *= frac;
        // record the fraction (out of 25) of shadow tests that are in shadow
        // and attenuate illumination accordingly

	    vec3 L = normalize(-spot_light_directions[i]);
		vec3 brdf_color = Phong_BRDF(L, V, N, diffuseColor, specularColor, specularExponent);
        //Lommel Seeligger Reflectance
        // vec3 brdf_color = LS_BRDF(L, V, N, diffuseColor);

        //Ward BRDF
        // vec3 brdf_color = Ward_BRDF(L, V, N, diffuseColor, specular_color);
        
        //Blinn-Phong Reflectance
        // vec3 H = normalize(L+V);
        // vec3 brdf_color = Blinn_Phong_BRDF(L, H, N, diffuseColor, specularColor, specularExponent);
	    
        //Oren-Nayar Reflectance
        // float roughness=0.5;
        // vec3 brdf_color = Oren_Nayar_BRDF(L, V, N, roughness, diffuseColor);
	    Lo += intensity * brdf_color;
    }

    fragColor = vec4(Lo, 1);
}



