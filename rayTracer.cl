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
    float roughness;
    float3 emission;
    float emissionIntensity;
    float clearCoat;
    float3 specular;
} Material;

Material createMaterial(float3 diffuse, float roughness, float clearCoat, float3 specular) {
    Material mat;
    mat.diffuse = diffuse;
    mat.roughness = roughness;
    mat.emission = (float3)(0.0f, 0.0f, 0.0f);
    mat.emissionIntensity = 0;
    mat.clearCoat = clearCoat;
    mat.specular = specular;
    return mat;
}

Material createEmissive(float3 emission, float emissionIntensity) {
    Material mat;
    mat.diffuse = (float3)(0.0f, 0.0f, 0.0f);
    mat.roughness = 0.0;
    mat.emission = emission;
    mat.emissionIntensity = emissionIntensity;
    mat.clearCoat = 0.0;
    mat.specular = (float3)(0.0f, 0.0f, 0.0f);
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

typedef struct {
    float3 pos;
    float3 scale;
    Material mat;
} Cube;

Cube createCube(float3 pos, float3 scale, Material mat) {
    Cube cube;
    cube.pos = pos;
    cube.scale = scale;
    cube.mat = mat;
    return cube;
}

float random(uint seed) {
    seed *= (seed + 195439) * (seed + 4877395) * (seed + 845921);
    return seed / 4294967296.0f;
}

float normalRandom(uint seed)
{
    float rand = random(seed);
    float theta = 2 * M_PI * rand;
    float rho = sqrt(-2 * log(rand));
    return rho * cos(theta);
}

float3 randomOnSphere(uint seed)
{
    return (float3)(normalRandom(seed + 893247), normalRandom(seed + 362793), normalRandom(seed + 9027126));
}

float2 randomOnCircle(uint seed)
{
    float angle = random(seed) * 2 * M_PI;
    float2 pointOnCircle = (float2)(cos(angle), sin(angle));
    return pointOnCircle * sqrt(random(seed + 398732));
}

float3 reflect(float3 dir, float3 normal)
{
    // r = d − 2(d⋅n)n
    return dir - normal * 2 * dot(dir, normal);
}

float3 slerp(float3 startDir, float3 endDir, float t) {
    // Normalize the start and end directions
    startDir = normalize(startDir);
    endDir = normalize(endDir);
    
    // Calculate the angle between the start and end directions
    float cosAngle = dot(startDir, endDir);
    float angle = acos(clamp(cosAngle, -1.0f, 1.0f)); // Clamp to avoid numerical issues
    
    // Handle the case when angle is small (linear interpolation fallback)
    if (fabs(angle) < 0.001f) {
        return mix(startDir, endDir, t); // Linear interpolation as fallback
    }

    // Calculate the slerp interpolation
    float sinAngle = sin(angle);
    float weightStart = sin((1 - t) * angle) / sinAngle;
    float weightEnd = sin(t * angle) / sinAngle;
    
    // Interpolate the directions using slerp
    return startDir * weightStart + endDir * weightEnd;
}

void intersectGroundPlane(Ray ray, RayHit *closestHit) {
    // Calculate intersection distance along the ray with the ground plane (y = 0)
    float t = -ray.pos.y / ray.dir.y;
    if (t > 0 && t < closestHit->dist) {
        closestHit->dist = t;
        closestHit->pos = ray.pos + ray.dir * t;
        closestHit->normal = (float3)(0.0f, 1.0f, 0.0f);
        // Create and assign a new material
        closestHit->mat = createMaterial((float3)(0.7f, 0.7f, 0.7f), 0.7, 0.2, (float3)(0.8f, 0.8f, 0.8f));
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

float3 cubeNormal(Ray ray, float3 hitPoint, Cube cube) {
    float3 p = cube.pos;
    float3 s = cube.scale;
    float epsilon = 0.0001f; // Adjust as needed for floating point precision

    if (fabs(hitPoint.x - (p.x - s.x / 2)) < epsilon) return (float3)(-1, 0, 0);
    else if (fabs(hitPoint.x - (p.x + s.x / 2)) < epsilon) return (float3)(1, 0, 0);
    else if (fabs(hitPoint.y - (p.y - s.y / 2)) < epsilon) return (float3)(0, -1, 0);
    else if (fabs(hitPoint.y - (p.y + s.y / 2)) < epsilon) return (float3)(0, 1, 0);
    else if (fabs(hitPoint.z - (p.z - s.z / 2)) < epsilon) return (float3)(0, 0, -1);
    else if (fabs(hitPoint.z - (p.z + s.z / 2)) < epsilon) return (float3)(0, 0, 1);

    return (float3)(0, 0, 0); // Should never happen; indicates an error
}

void intersectCube(Ray ray, RayHit *closestHit, Cube cube) {
    float3 p = cube.pos - cube.scale / 2.0f;
    float3 s = cube.scale;

    // Calculate the minimum and maximum t values for each axis
    float tminX = (p.x - ray.pos.x) / ray.dir.x;
    float tmaxX = (p.x + s.x - ray.pos.x) / ray.dir.x;
    float tminY = (p.y - ray.pos.y) / ray.dir.y;
    float tmaxY = (p.y + s.y - ray.pos.y) / ray.dir.y;
    float tminZ = (p.z - ray.pos.z) / ray.dir.z;
    float tmaxZ = (p.z + s.z - ray.pos.z) / ray.dir.z;

    // Calculate the actual minimum and maximum t values
    float tmin = fmax(fmax(fmin(tminX, tmaxX), fmin(tminY, tmaxY)), fmin(tminZ, tmaxZ));
    float tmax = fmin(fmin(fmax(tminX, tmaxX), fmax(tminY, tmaxY)), fmax(tminZ, tmaxZ));

    // Check if the ray misses the cube
    if (tmax < 0 || tmin > tmax || tmin > closestHit->dist) return;

    // Check if the ray hits the cube
    if (tmin > 0 && tmin < closestHit->dist) {
        closestHit->dist = tmin;
        closestHit->pos = ray.pos + ray.dir * tmin;
        closestHit->normal = cubeNormal(ray, closestHit->pos, cube); // Assuming cubeNormal function exists and is adapted below
        closestHit->mat = cube.mat; // Direct assignment assuming Material can be copied
    }
}

RayHit trace(Ray ray) {
    RayHit closestHit;
    closestHit.dist = FLT_MAX;
    
    intersectGroundPlane(ray, &closestHit);

    float3 color1 = (float3)(121, 224, 238) / 255;
    float3 color2 = (float3)(152, 238, 204) / 255;
    float3 color3 = (float3)(208, 245, 190) / 255;
    float3 color4 = (float3)(251, 255, 220) / 255;
    
    intersectSphere(ray, &closestHit, createSphere((float3)(0, 3, 7), 3, createMaterial(color1, 1, 0, (float3)(1, 1, 1))));
    intersectSphere(ray, &closestHit, createSphere((float3)(-2.5, 1, 3.5), 1, createMaterial(color2, 1, 0, (float3)(1, 1, 1))));
    intersectSphere(ray, &closestHit, createSphere((float3)(3, 1.5, 3), 1.5, createMaterial(color3, 1, 0, (float3)(1, 1, 1))));
    intersectSphere(ray, &closestHit, createSphere((float3)(10, 4, 10), 4, createMaterial(color4, 1, 0, (float3)(1, 1, 1))));
    intersectSphere(ray, &closestHit, createSphere((float3)(-7, 2, 6), 2, createMaterial(color1, 1, 0, (float3)(1, 1, 1))));
    intersectSphere(ray, &closestHit, createSphere((float3)(-6.5, 1.5, 12), 1.5, createMaterial(color2, 1, 0, (float3)(1, 1, 1))));
    
    return closestHit;
}

void shineRay(Ray *ray, RayHit hit, uint seed) {
    // Handle light emission
    if (length(hit.mat.emission) > 0) {
        ray->carriedLight = hit.mat.emission * ray->energy * hit.mat.emissionIntensity;
        ray->energy = (float3)(0.0f, 0.0f, 0.0f);
        return;
    }

    if (hit.dist < FLT_MAX) {
        // Offset to correct for precision issues
        ray->pos = hit.pos + hit.normal * 0.001f;
        
        // Specularity and reflection calculations
        bool isSpecular = (random(seed) <= hit.mat.clearCoat);
        
        float3 reflectedRay = reflect(ray->dir, hit.normal);
        float3 randomRay = normalize(hit.normal + randomOnSphere(seed));
        float3 diffuseRay = slerp(reflectedRay, randomRay, hit.mat.roughness);

        ray->dir = reflectedRay;
        ray->dir = isSpecular ? reflectedRay : diffuseRay;
        float3 reflectColor = isSpecular ? hit.mat.specular : hit.mat.diffuse;
        
        ray->energy *= reflectColor;

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

Ray createCameraRay(float3 cameraPos, float2 uv, float fov, float aperture, float focalLength, uint seed) {
    Ray ray;
    ray.energy = (float3)(1.0f, 1.0f, 1.0f);

    // Generate aperture offset
    float2 rand = randomOnCircle(seed);
    float3 rayPos = (float3)(rand.x * aperture, rand.y * aperture, 0.0f);

    // Set ray position and direction
    ray.pos = cameraPos + rayPos;
    float3 dirComponent = (float3)(uv.x * fov * focalLength, uv.y * fov * focalLength, focalLength);
    ray.dir = normalize(dirComponent - rayPos);

    return ray;
}

inline int float3_color(float3 color) {
    // Convert each component to an integer in the range [0, 255]
    color = clamp(color, 0.0f, 1.0f);

    int ir = (int)(color.x * 255.0f);
    int ig = (int)(color.y * 255.0f);
    int ib = (int)(color.z * 255.0f);

    // Combine the components into a single int color in ARGB format
    int colorInt = (255 << 24) | (ir << 16) | (ig << 8) | ib;

    return colorInt;
}

__kernel void rayTracer(__global int* pixelData, const int width, const int height) {
    int reflectionLimit = 1000;
    int samples = 1000;
    
    int i = get_global_id(0);
    int j = get_global_id(1);

    float u = (float)i / height * 2.0f - 1.0f;
    float v = (float)j / width * 2.0f - 1.0f;
    u *= width / height;
    v = -v;

    static float3 camera = (float3)(0.0f, 3.0f, 0.0f);
    static float fov = 1;
    static float aperture = 0;
    static float focalLength = 6;

    float3 light = (float3)(0.0f, 0.0f, 0.0f);

    for (int s = 0; s < samples; s++)
    {
        uint sampleId = s + samples * (j + width * i);
        Ray ray = createCameraRay(camera,  (float2)(u, v), fov, aperture, focalLength, sampleId);
        for (int r = 0; r < reflectionLimit; r++)
        {
            RayHit hit = trace(ray);
            shineRay(&ray, hit, sampleId);
            if (length(ray.energy) <= 0.001) break;                   
        }

        float weight = 1 / (float)(s + 1);
        light = mix(light, ray.carriedLight, weight);
    }

    int pixel = float3_color(light);

    pixelData[j * width + i] = pixel;
}