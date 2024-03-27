import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import com.aparapi.Kernel;

public class RayTracer extends RayMath
{   
    public class RayTracingKernel extends Kernel {
        private final int[] pixelArray;
        private final int width;
        private final int height;

        public RayTracingKernel(int[] pixelArray, int width, int height) {
            this.pixelArray = pixelArray;
            this.width = width;
            this.height = height;
        }

        @Override
        public void run() {
            int gid = getGlobalId();
            int x = gid % width;
            int y = gid / width;

            // Cast a ray into the scene and calculate the color of the pixel
            // This is where you would put your ray tracing code
            // pixelArray[gid] = castRay(x, y);
        }
    }

    static int width = 1000;
    static int height = 1000;
    
    static int reflectionLimit = 1000;
    static int samples = 4000;
    
    static Vector3 camera = new Vector3(0, 3, 0);
    static double fov = 1;
    static double aperture = 0.6;
    static double focalLength = 6;

    public static void main(String[] args) throws IOException
    {
        System.out.println("Rendering... (" + width + "px, " + height + "px, " + samples + " samples)");
        Vector3[] pixelArray = new Vector3[width * height];
        for (int i = 0; i < pixelArray.length; i++)
        {
            pixelArray[i] = new Vector3();
        }

        for (int s = 0; s < samples; s++)
        {            
            BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
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
                    Vector3 average = pixelArray[i * height + j].mul(1 - weight).add(ray.carriedLight.mul(weight));
                    pixelArray[i * height + j] = average;
                    average = average.mul(255);
                    average = new Vector3(Clamp(average.x, 0, 255), Clamp(average.y, 0, 255), Clamp(average.z, 0, 255));
                    int pixel = Color(average);
                    img.setRGB(j, i, pixel);
                }
            }

            // Save as PNG
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