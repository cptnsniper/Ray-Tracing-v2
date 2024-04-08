typedef struct {
    float3 pos;
    float3 dir;
    float3 energy;
    float3 carriedLight;
} Ray;

Ray createRay(float3 pos, float3 dir, float3 energy, float3 carriedLight) {
    Ray ray;
    ray.pos = pos;
    ray.dir = dir;
    ray.energy = energy;
    ray.carriedLight = carriedLight;
    return ray;
}

typedef struct {
    float3 diffuse;
    // double roughness;
    // // float3 emission;
    // // double emissionIntensity;
    // // double clearCoat;
    // // float3 specular;
} Material;

Material createMaterial(float3 diffuse) {
    Material mat;
    mat.diffuse = diffuse;
    // mat.roughness = roughness;
    // mat.emission = emission;
    // mat.emissionIntensity = emissionIntensity;
    // mat.clearCoat = clearCoat;
    // mat.specular = specular;
    return mat;
}

typedef struct {
    float3 pos;
    float dist;
    float3 normal;
    Material mat;
} RayHit;

RayHit createRayHit(float3 pos, float dist, float3 normal, Material mat) {
    RayHit rayHit;
    rayHit.pos = pos;
    rayHit.dist = dist;
    rayHit.normal = normal;
    rayHit.mat = mat;
    return rayHit;
}

typedef struct {
    float3 pos;
    float radius;
    Material mat;
} Sphere;

Sphere createSphere(float3 pos, float radius, Material mat) {
    Sphere sphere;
    sphere.pos = pos;
    sphere.radius = radius;
    sphere.mat = mat;
    return sphere;
}

void IntersectGroundPlane(Ray ray, RayHit *closestHit) {
    // Calculate intersection distance along the ray with the ground plane (y = 0)
    float t = -ray.pos.y / ray.dir.y;
    if (t > 0 && t < closestHit->dist) {
        closestHit->dist = t;
        closestHit->pos = ray.pos + ray.dir * t;
        closestHit->normal = (float3)(0.0f, 1.0f, 0.0f);
        // Create and assign a new material
        closestHit->mat = createMaterial((float3)(0.0f, 1.0f, 0.0f));
    }
}

void intersectSphere(Ray ray, RayHit *closestHit, Sphere sphere) {
    float3 d = ray.pos - sphere.pos;
    float p1 = -dot(ray.dir, d);
    float p2sqr = pow(p1, 2) - dot(d, d) + pow(sphere.radius, 2);
    if (p2sqr < 0) return;
    float p2 = sqrt(p2sqr);
    float t = p1 - p2 > 0 ? p1 - p2 : p1 + p2;
    if (t > 0 && t < closestHit->dist) {
        closestHit->dist = t;
        closestHit->pos = ray.pos + ray.dir * t;
        closestHit->normal = normalize(closestHit->pos - sphere.pos);
        closestHit->mat = sphere.mat;
    }
}

RayHit trace(Ray ray) {
    RayHit closestHit = createRayHit((float3)(0.0f, 0.0f, 0.0f), FLT_MAX, (float3)(0.0f, 0.0f, 0.0f), createMaterial((float3)(0.0f, 0.0f, 0.0f)));
    
    IntersectGroundPlane(ray, &closestHit);
    intersectSphere(ray, &closestHit, createSphere((float3)(0.0f, 4.0f, 5.0f), 0.5f, createMaterial((float3)(1.0f, 0.0f, 0.0f))));
    
    return closestHit;
}

void shineRay(Ray *ray, RayHit hit) {
    // Handle light emission
    // if (length(hit.mat.emission) > 0) {
    //     ray.carriedLight = hit.mat.emission * hit.mat.emissionIntensity * ray.energy;
    //     ray.energy = (float3)(0.0f, 0.0f, 0.0f);
    //     return ray;
    // }

    if (hit.dist < FLT_MAX) {
        // Offset to correct for precision issues
        // ray.pos = hit.pos + hit.normal * 0.001f;
        
        // // Specularity and reflection calculations
        // bool isSpecular = (frand() <= hit.mat.clearCoat); // frand() needs to be defined or replaced with an appropriate RNG function
        
        // float3 reflectedRay = Reflect(ray.dir, hit.normal);
        // float3 randomRay = normalize(hit.normal + RandomOnSphere());
        // float3 diffuseRay = Slerp(reflectedRay, randomRay, hit.mat.roughness);

        // ray.dir = isSpecular ? reflectedRay : diffuseRay;
        // float3 reflectColor = isSpecular ? hit.mat.specular : hit.mat.diffuse;
        
        // ray.energy = ray.energy * reflectColor;

        ray->energy = (float3)(0.0f, 0.0f, 0.0f);
        ray->carriedLight = hit.mat.diffuse;

    } else {
        // Sky and sun handling
        float3 sunDirection = normalize((float3)(-1.0f, 1.0f, 1.0f));
        if (length(ray->dir - sunDirection) < 0.2f) {
            ray->carriedLight = (float3)(1.0f, 0.84f, 0.67f) * 20.0f;
            ray->energy = (float3)(0.0f, 0.0f, 0.0f);
            return;
        }

        float weight = ray->dir.y / 2.0f + 0.5f;
        float3 bottomColor = (float3)(1.0f, 1.0f, 1.0f);
        float3 topColor = (float3)(0.83f, 0.92f, 1.0f);
        ray->carriedLight = (bottomColor * (1.0f - weight) + topColor * weight) * ray->energy * 0.9f;
        ray->energy = (float3)(0.0f, 0.0f, 0.0f);
    }
}

Ray createCameraRay(float3 cameraPos, float2 uv, float fov, float aperture, float focalLength) {
    Ray ray;

    // Generate aperture offset
    float2 randomOnCircle = (float2)(0.0f, 0.0f); // RandomOnCircle() needs to be defined or replaced with an appropriate RNG function
    float3 rayPos = (float3)(randomOnCircle.x, randomOnCircle.y, 0.0f);

    // Set ray position and direction
    ray.pos = cameraPos + rayPos; // Vector addition
    float3 dirComponent = (float3)(uv.x * fov * focalLength, uv.y * fov * focalLength, focalLength);
    ray.dir = normalize(dirComponent - rayPos); // Vector subtraction and normalization

    return ray;
}

inline int float3_color(float3 color) {
    // Convert each component to an integer in the range [0, 255]
    int ir = (int)(color.x * 255.0f);
    int ig = (int)(color.y * 255.0f);
    int ib = (int)(color.z * 255.0f);
    // Combine the components into a single int color in ARGB format
    int colorInt = (255 << 24) | (ir << 16) | (ig << 8) | ib;

    return colorInt;
}

__kernel void rayTracer(__global int* pixelData, const int width, const int height) {
    static int reflectionLimit = 1000;
    
    int i = get_global_id(0);
    int j = get_global_id(1);

    float u = (float)i / height * 2.0f - 1.0f;
    float v = (float)j / width * 2.0f - 1.0f;
    u *= width / height;
    v = 1.0f - v;

    static float3 camera = (float3)(0.0f, 3.0f, 0.0f);
    static double fov = 1;
    static double aperture = 0.6;
    static double focalLength = 1;

    // Ray ray = createCameraRay(camera,  (float2)(u, v), fov, aperture, focalLength);
    Ray ray = createRay(camera, (float3)(u, v, 1.0f), (float3)(1.0f, 1.0f, 1.0f), (float3)(0.0f, 0.0f, 0.0f));
    // for (int r = 0; r < reflectionLimit; r++)
    // {
        RayHit hit = trace(ray);
        shineRay(&ray, hit);
    //     if (length(ray.energy) <= 0.001) break;                   
    // }
    
    int pixel = float3_color(ray.carriedLight);
    pixelData[j * width + i] = float3_color(pixel);
}