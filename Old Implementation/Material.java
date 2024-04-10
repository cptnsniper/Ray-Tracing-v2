public class Material
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