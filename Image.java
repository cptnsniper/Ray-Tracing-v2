import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;

public class Image
{
    public static class Vector3
    {
        private double x;
        private double y;
        private double z;
        public Vector3()
        {
            this.x = 0;
            this.y = 0;
            this.z = 0;
        } 
        public Vector3(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }
        @Override
        public String toString()
        {
            return "(" + x + ", " + y + ", " + z + ")";
        }
        
        public static double dist(Vector3 vect1, Vector3 vect2)
        {
            // D = sqrt(Δx^2 + Δy^2 + Δz^2)
            return Math.sqrt(Math.pow(vect2.x - vect1.x, 2) + Math.pow(vect2.y - vect1.y, 2) + Math.pow(vect2.z - vect1.z, 2));
        }
        public double mag()
        {
            return dist(new Vector3(), this);
        }
        public Vector3 normal()
        {
            double mag = this.mag();
            return new Vector3(this.x / mag, this.y / mag, this.z / mag);
        }
        public static double dot(Vector3 vect1, Vector3 vect2)
        {
            return vect1.x * vect2.x + vect1.y * vect2.y + vect1.z * vect2.z;
        }
        public static Vector3 mul(Vector3 vect1, Vector3 vect2)
        {
            return new Vector3(vect1.x * vect2.x, vect1.y * vect2.y, vect1.z * vect2.z);
        }
        public static Vector3 div(Vector3 vect1, Vector3 vect2)
        {
            return new Vector3(vect1.x / vect2.x, vect1.y / vect2.y, vect1.z / vect2.z);
        }
        public static Vector3 slerp(Vector3 startDir, Vector3 endDir, double t)
        {
            startDir = startDir.normal();
            endDir = endDir.normal();
        
            // Calculate the angle between the start and end directions
            double angle = Math.acos(Vector3.dot(startDir, endDir));
        
            // Calculate the slerp interpolation
            double sinAngle = Math.sin(angle);
            double weightStart = Math.sin((1 - t) * angle) / sinAngle;
            double weightEnd = Math.sin(t * angle) / sinAngle;
        
            // Interpolate the directions using slerp
            Vector3 lerpedDir = startDir.mul(weightStart).add(endDir.mul(weightEnd));
        
            lerpedDir = lerpedDir.normal();
        
            return lerpedDir;
        }
        
        public Vector3 mul(double factor)
        {
            return new Vector3(this.x * factor, this.y * factor, this.z * factor);
        }
        public Vector3 div(double factor)
        {
            return new Vector3(this.x / factor, this.y / factor, this.z / factor);
        }
        public Vector3 add(Vector3 vector)
        {
            return new Vector3(this.x + vector.x, this.y + vector.y, this.z + vector.z);
        }
        public Vector3 sub(Vector3 vector)
        {
            return new Vector3(this.x - vector.x, this.y - vector.y, this.z - vector.z);
        }
    }
    
    public static class Vector2
    {
        private double x;
        private double y;
        public Vector2()
        {
            this.x = 0;
            this.y = 0;
        } 
        public Vector2(double x, double y)
        {
            this.x = x;
            this.y = y;
        }
        @Override
        public String toString()
        {
            return "(" + x + ", " + y + ")";
        }
        
        public static double dist(Vector2 vect1, Vector2 vect2)
        {
            // D = sqrt(Δx^2 + Δy^2)
            return Math.sqrt(Math.pow(vect2.x - vect1.x, 2) + Math.pow(vect2.y - vect1.y, 2));
        }
        public double mag()
        {
            return dist(new Vector2(), this);
        }
        public Vector2 normal()
        {
            double mag = this.mag();
            return new Vector2(this.x / mag, this.y / mag);
        }
        public static double dot(Vector2 vect1, Vector2 vect2)
        {
            return vect1.x * vect2.x + vect1.y * vect2.y;
        }
        public static Vector2 mul(Vector2 vect1, Vector2 vect2)
        {
            return new Vector2(vect1.x * vect2.x, vect1.y * vect2.y);
        }
        
        public Vector2 mul(double factor)
        {
            return new Vector2(this.x * factor, this.y * factor);
        }
        public Vector2 div(double factor)
        {
            return new Vector2(this.x / factor, this.y / factor);
        }
        public Vector2 add(Vector2 vector)
        {
            return new Vector2(this.x + vector.x, this.y + vector.y);
        }
        public Vector2 sub(Vector2 vector)
        {
            return new Vector2(this.x - vector.x, this.y - vector.y);
        }
    }
    
    public static class Ray
    {
        Vector3 pos;
        Vector3 dir;
        Vector3 energy;
        Vector3 carriedLight;
        public Ray()
        {
            this.pos = new Vector3();
            this.dir = new Vector3();
            this.energy = new Vector3(1, 1, 1);
            this.carriedLight = new Vector3();
        }
        public Ray(Vector3 pos, Vector3 dir)
        {
            this.pos = pos;
            this.dir = dir;
            this.energy = new Vector3(1, 1, 1);
            this.carriedLight = new Vector3();
        }
        
        public Vector3 reflect(Vector3 dir, Vector3 normal)
        {
            // r = d − 2(d⋅n)n
            return dir.sub(normal.mul(2 * Vector3.dot(dir, normal)));
        }
    }
    
    public static class RayHit
    {
        Vector3 pos;
        double dist;
        Vector3 normal;
        Material mat;
        public RayHit()
        {
            this.pos = new Vector3();
            this.dist = Double.POSITIVE_INFINITY;
            this.normal = new Vector3();
            this.mat = new Material(new Vector3(), 0);
        }
        public RayHit(Vector3 pos, double dist, Vector3 normal, Material mat)
        {
            this.pos = pos;
            this.dist = dist;
            this.normal = normal;
            this.mat = mat;
        }
    }
    
    public static class Material
    {
        Vector3 diffuse;
        double roughness;
        Vector3 emission;
        double emissionIntensity;
        double clearCoat;
        Vector3 specular;
        public Material(Vector3 diffuse, double roughness, double clearCoat, Vector3 specular)
        {
            this.diffuse = diffuse;
            this.roughness = roughness;
            this.emission = new Vector3();
            this.emissionIntensity = 0;
            this.clearCoat = clearCoat;
            this.specular = specular;
        }
        public Material(Vector3 emission, double intensity)
        {
            this.emission = emission;
            this.emissionIntensity = intensity;
        }
    }
    
    public static class Sphere
    {
        Vector3 pos;
        double radius;
        Material mat;
        public Sphere(Vector3 pos, double radius, Material mat)
        {
            this.pos = pos;
            this.radius = radius;
            this.mat = mat;
        }
    }
    public static class Cube
    {
        Vector3 pos;
        Vector3 scale;
        Material mat;
        public Cube(Vector3 pos, Vector3 scale, Material mat)
        {
            this.pos = pos;
            this.scale = scale;
            this.mat = mat;
        }
    }
    
    public static RayHit IntersectGroundPlane(Ray ray, RayHit closestHit)
    {
        // Calculate distance along the ray where the ground plane is intersected
        double t = -ray.pos.y / ray.dir.y;
        if (t > 0 && t < closestHit.dist)
        {
            closestHit.dist = t;
            closestHit.pos = ray.pos.add(ray.dir.mul(t));
            closestHit.normal = new Vector3(0.0f, 1.0f, 0.0f);
            
            closestHit.mat = new Material(new Vector3(0.8, 0.8, 0.8), 1, 0, new Vector3(1, 1, 1));
            
            // Checkerboard
            // Material white = new Material(new Vector3(0.8, 0.8, 0.8), 0.7, 0.2, new Vector3(1, 1, 1));
            // Material black = new Material(new Vector3(0.2, 0.2, 0.2), 0.7, 0.2, new Vector3(1, 1, 1));
            // if (Math.round(Math.abs(closestHit.pos.x)) % 2 == 0)
            // {
                // if (Math.round(Math.abs(closestHit.pos.z)) % 2 == 0)
                // {
                    // closestHit.mat = black;
                // }
                // else
                // {
                    // closestHit.mat = white;
                // }
            // }
            // else
            // {
                // if (Math.round(Math.abs(closestHit.pos.z)) % 2 == 1)
                // {
                    // closestHit.mat = black;
                // }
                // else
                // {
                    // closestHit.mat = white;
                // }
            // }
            return closestHit;
        }
        return closestHit;
    }
    
    public static RayHit IntersectSphere(Ray ray, RayHit closestHit, Sphere sphere)
    {
        // Calculate distance along the ray where the sphere is intersected
        Vector3 d = ray.pos.sub(sphere.pos);
        double p1 = -Vector3.dot(ray.dir, d);
        double p2sqr = Math.pow(p1, 2) - Vector3.dot(d, d) + Math.pow(sphere.radius, 2);
        if (p2sqr < 0) return closestHit;
        double p2 = Math.sqrt(p2sqr);
        double t = p1 - p2 > 0 ? p1 - p2 : p1 + p2;
        if (t > 0 && t < closestHit.dist)
        {
            closestHit.dist = t;
            closestHit.pos = ray.pos.add(ray.dir.mul(t));
            closestHit.normal = closestHit.pos.sub(sphere.pos).normal();
            closestHit.mat = sphere.mat;
            return closestHit;
        }
        return closestHit;
    }
    
    public static RayHit IntersectCube(Ray ray, RayHit closestHit, Cube cube)
    {
        Vector3 p = cube.pos.sub(cube.scale.div(2));
        Vector3 s = cube.scale;
    
        // Calculate the minimum and maximum t values for each axis
        double tminX = (p.x - ray.pos.x) / ray.dir.x;
        double tmaxX = (p.x + s.x - ray.pos.x) / ray.dir.x;
        double tminY = (p.y - ray.pos.y) / ray.dir.y;
        double tmaxY = (p.y + s.y - ray.pos.y) / ray.dir.y;
        double tminZ = (p.z - ray.pos.z) / ray.dir.z;
        double tmaxZ = (p.z + s.z - ray.pos.z) / ray.dir.z;
    
        // Calculate the actual minimum and maximum t values
        double tmin = Math.max(Math.max(Math.min(tminX, tmaxX), Math.min(tminY, tmaxY)), Math.min(tminZ, tmaxZ));
        double tmax = Math.min(Math.min(Math.max(tminX, tmaxX), Math.max(tminY, tmaxY)), Math.max(tminZ, tmaxZ));
    
        // Check if the ray misses the cube
        if (tmax < 0 || tmin > tmax || tmin > closestHit.dist) return closestHit;
    
        // Check if the ray hits the cube
        if (tmin > 0 && tmin < closestHit.dist)
        {
            closestHit.dist = tmin;
            closestHit.pos = ray.pos.add(ray.dir.mul(tmin));
            closestHit.normal = cubeNormal(ray, closestHit.pos, cube);
            closestHit.mat = cube.mat;
            // closestHit.mat = new Material(new Vector3(1, 1, 1), 1, new Vector3(1, 1, 1), 1);
            return closestHit;
        }
    
        return closestHit;
    }
    
    public static Vector3 cubeNormal(Ray ray, Vector3 hitPoint, Cube cube)
    {
        Vector3 p = cube.pos;
        Vector3 s = cube.scale;
    
        double epsilon = 0.0001;
    
        if (Math.abs(hitPoint.x - (p.x - s.x / 2)) < epsilon)
        {
            return new Vector3(-1, 0, 0);
        }
        else if (Math.abs(hitPoint.x - (p.x + s.x / 2)) < epsilon)
        {
            return new Vector3(1, 0, 0);
        }
        else if (Math.abs(hitPoint.y - (p.y - s.y / 2)) < epsilon)
        {
            return new Vector3(0, -1, 0);
        }
        else if (Math.abs(hitPoint.y - (p.y + s.y / 2)) < epsilon)
        {
            return new Vector3(0, 1, 0);
        }
        else if (Math.abs(hitPoint.z - (p.z - s.z / 2)) < epsilon)
        {
            return new Vector3(0, 0, -1);
        }
        else if (Math.abs(hitPoint.z - (p.z + s.z / 2)) < epsilon)
        {
            return new Vector3(0, 0, 1);
        }
    
        return new Vector3();
    }
    
    public static RayHit Trace(Ray ray)
    {
        RayHit closestHit = new RayHit();
        closestHit = IntersectGroundPlane(ray, closestHit);
        // Corenell Box
        // Material whiteWall = new Material(new Vector3(0.8, 0.8, 0.8), 0.7, 0, new Vector3(1, 1, 1));
        // Material redWall = new Material(new Vector3(0.8, 0.3, 0.3), 0.7, 0, new Vector3(1, 1, 1));
        // Material blueWall = new Material(new Vector3(0.3, 0.3, 0.8), 0.7, 0, new Vector3(1, 1, 1));
        
        // closestHit = IntersectCube(ray, closestHit, new Cube(new Vector3(0, 10, 10), new Vector3(22, 0.1, 22), whiteWall));
        // closestHit = IntersectCube(ray, closestHit, new Cube(new Vector3(0, 10, 20), new Vector3(22, 22, 0.1), whiteWall));
        // closestHit = IntersectCube(ray, closestHit, new Cube(new Vector3(10, 10, 10), new Vector3(0.1, 22, 22), blueWall));
        // closestHit = IntersectCube(ray, closestHit, new Cube(new Vector3(-10, 10, 10), new Vector3(0.1, 22, 22), redWall));
        // closestHit = IntersectCube(ray, closestHit, new Cube(new Vector3(0, 10, -1), new Vector3(22, 22, 0.1), whiteWall));
        // Light
        // closestHit = IntersectCube(ray, closestHit, new Cube(new Vector3(0, 9.9, 10), new Vector3(8, 0.1, 8), new Material(new Vector3(1, 0.82, 0.64), 10)));
        
        Vector3 color1 = new Vector3(121, 224, 238).div(255);
        Vector3 color2 = new Vector3(152, 238, 204).div(255);
        Vector3 color3 = new Vector3(208, 245, 190).div(255);
        Vector3 color4 = new Vector3(251, 255, 220).div(255);
        
        closestHit = IntersectSphere(ray, closestHit, new Sphere(new Vector3(0, 3, 7), 3, new Material(color1, 1, 0.2, new Vector3(1, 1, 1))));
        closestHit = IntersectSphere(ray, closestHit, new Sphere(new Vector3(-2, 0.5, 4), 1, new Material(color2, 1, 0.2, new Vector3(1, 1, 1))));
        closestHit = IntersectSphere(ray, closestHit, new Sphere(new Vector3(3, 1.5, 3), 1.5, new Material(color3, 1, 0.2, new Vector3(1, 1, 1))));
        closestHit = IntersectSphere(ray, closestHit, new Sphere(new Vector3(8, 4, 8), 4, new Material(color4, 1, 0.2, new Vector3(1, 1, 1))));
        closestHit = IntersectSphere(ray, closestHit, new Sphere(new Vector3(-7, 2, 6), 2, new Material(color1, 1, 0.2, new Vector3(1, 1, 1))));
        closestHit = IntersectSphere(ray, closestHit, new Sphere(new Vector3(-6.5, 1.5, 12), 1.5, new Material(color2, 1, 0.2, new Vector3(1, 1, 1))));
        
        // closestHit = IntersectSphere(ray, closestHit, new Sphere(new Vector3(-20, 40, 100), 25, new Material(new Vector3(1, 1, 1), 1, new Vector3(1, 1, 0.8), 7)));
        
        return closestHit;
    }
    
    public static Ray ShineRay(Ray ray, RayHit hit)
    {
        // If the ray hits a light
        if (hit.mat.emission.mag() > 0)
        {
            ray.carriedLight = Vector3.mul(hit.mat.emission.mul(hit.mat.emissionIntensity), ray.energy);
            ray.energy = new Vector3();
            return ray;
        }
        // If the ray does not fly into the sky
        if (hit.dist < Double.POSITIVE_INFINITY)
        {   
            // Offset to correct for double imprecision
            ray.pos = hit.pos.add(hit.normal.mul(0.001));
            
            boolean isSpecular = Math.random() <= hit.mat.clearCoat;
            
            // Calculating reflection diriction
            Vector3 reflectedRay = ray.reflect(ray.dir, hit.normal);
            Vector3 randomRay = hit.normal.add(RandomOnSphere()).normal();
            Vector3 diffuseRay = Vector3.slerp(reflectedRay, randomRay, hit.mat.roughness);
            
            // If specular use specular dir and color else use the diffuse ones
            ray.dir = isSpecular ? reflectedRay : diffuseRay;
            Vector3 reflectColor = isSpecular ? hit.mat.specular : hit.mat.diffuse;
            
            ray.energy = Vector3.mul(ray.energy, reflectColor);
            
            // Return Normal
            // ray.energy = new Vector3();
            // ray.carriedLight = randomRay.div(2).add(new Vector3(0.5, 0.5, 0.5)).mul(255);
            
            return ray;
        }
        else
        {
            // Calculate if the ray hits the sun
            Vector3 sun = new Vector3(-1, 1, 1).normal();
            if (ray.dir.sub(sun).mag() < 0.2)
            {
                ray.carriedLight = new Vector3(1, 0.84, 0.67).mul(20);
                ray.energy = new Vector3();
                return ray;
            }
            
            // Lerp between top and bottem colors for the sky
            double weight = ray.dir.normal().y / 2 + 0.5;
            Vector3 bottemColor = new Vector3(1, 1, 1);
            Vector3 topColor = new Vector3(0.83, 0.92, 1);
            ray.carriedLight = Vector3.mul(bottemColor.mul(1 - weight).add(topColor.mul(weight)), ray.energy.mul(0.9));
            ray.energy = new Vector3();
            
            return ray;
        }
    }
    
    public static Ray CreateCameraRay(Vector3 cameraPos, Vector2 uv, double fov, double aperture, double focalLength)
    {
        Ray ray = new Ray();
        // Aperture offset
        Vector2 randomOnCircle = RandomOnCircle();
        Vector3 rayPos = new Vector3(randomOnCircle.x * aperture, randomOnCircle.y * aperture, 0);
        ray.pos = cameraPos.add(rayPos);
        ray.dir = new Vector3(uv.x * fov * focalLength, uv.y * fov * focalLength, focalLength).sub(rayPos).normal();
        return ray;
    }

    public static int Color(Vector3 rgb)
    {
        rgb.x = ((int)Math.rint(rgb.x) << 16) & 0x00FF0000;
        rgb.y = ((int)Math.rint(rgb.y) << 8) & 0x0000FF00;
        rgb.z = (int)Math.rint(rgb.z) & 0x000000FF;
    
        return 0xFF000000 | (int)rgb.x | (int)rgb.y | (int)rgb.z;
    }
    
    public static double Map(double value, double inMin, double inMax, double outMin, double outMax)
    {
        return outMin + ((outMax - outMin) / (inMax - inMin)) * (value - inMin);
    }
    
    public static double Clamp(double value, double min, double max)
    {
        return Math.min(Math.max(min, value), max);
    }
    
    public static double NormalRandom()
    {
        double random = Math.random();
        double theta = 2 * Math.PI * random;
        double rho = Math.sqrt(-2 * Math.log(random));
        return rho * Math.cos(theta);
    }
    public static Vector3 RandomOnSphere()
    {
        return new Vector3(NormalRandom(), NormalRandom(), NormalRandom());
    }
    public static Vector2 RandomOnCircle()
    {
        double angle = Math.random() * 2 * Math.PI;
        Vector2 pointOnCircle = new Vector2(Math.cos(angle), Math.sin(angle));
        return pointOnCircle.mul(Math.sqrt(Math.random()));
    }
    
    public static void main(String[] args) throws IOException
    {
        int width = 1000;
        int height = 1000;
        
        int reflectionLimit = 1000;
        int samples = 4000;
        
        Vector3 camera = new Vector3(0, 3, 0);
        double fov = 1;
        double aperture = 0.6;
        double focalLength = 6;
 
        System.out.println("Rendering... (" + width + "px, " + height + "px, " + samples + " samples)");
        Vector3[] pixelArray = new Vector3[width * height];
        for (int i = 0; i < pixelArray.length; i++)
        {
            pixelArray[i] = new Vector3();
        }

        for (int s = 0; s < samples; s++)
        {            
            BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            // BufferedImage toneMappedImg = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    double u = (double)j / width * 2 - 1;
                    double v = (double)i / height * 2 - 1;
                    u *= width / height;
                    v = -v;

                    Ray ray = CreateCameraRay(camera, new Vector2(u, v), fov, aperture, focalLength);
                    for (int r = 0; r < reflectionLimit; r++)
                    {
                        RayHit hit = Trace(ray);
                        ray = ShineRay(ray, hit);
                        if (ray.energy.mag() <= 0.001) break;                   
                    }
                    
                    double weight = 1 / (double)(s + 1);
                    // Color oldPixel = new Color(oldImg.getRGB(j, i));
                    // Vector3 oldResult = new Vector3(oldPixel.getRed(), oldPixel.getGreen(), oldPixel.getBlue());
                    Vector3 average = pixelArray[i * height + j].mul(1 - weight).add(ray.carriedLight.mul(weight));
                    pixelArray[i * height + j] = average;
                    average = average.mul(255);
                    average = new Vector3(Clamp(average.x, 0, 255), Clamp(average.y, 0, 255), Clamp(average.z, 0, 255));
                    Vector3 color = average;
                    int pixel = Color(average);
                    img.setRGB(j, i, pixel);
                    
                    // Tone Mapping
                    // color = color.div(4);
                    // double white = 2;
                    // double exposure = 10;
                    // color = color.mul(white * exposure);
                    // color = Vector3.div(Vector3.mul(color, new Vector3(1, 1, 1).add(color.div(white / white))), (new Vector3(1, 1, 1).add(color)));
                    // color =  new Vector3(Clamp(color.x, 0, 255), Clamp(color.y, 0, 255), Clamp(color.z, 0, 255));
                    // pixel = Color(color);
                    // toneMappedImg.setRGB(j, i, pixel);
                }
            }

            // Save as PNG
            // File file = new File("image\\image.png");
            // ImageIO.write(toneMappedImg, "png", file); 
            File file = new File("image\\Render Result [" + s + "].png");
            ImageIO.write(img, "png", file);
            
            int barLength = 20;
            double completePercent = (s + 1) / (double)samples;
            String bar = "";
            for (int b = 0; b < Math.round(completePercent * barLength); b++)
            {
                bar += "#";
            }
            for (int b = 0; b < barLength - Math.round(completePercent * barLength); b++)
            {
                bar += "-";
            }
            System.out.println("|" + bar + "| " + "Sample: (" + (s + 1) + " / " + samples + ")");

        }
        System.out.println("Rendering Finnished");
    }
}