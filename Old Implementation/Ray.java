public class Ray
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