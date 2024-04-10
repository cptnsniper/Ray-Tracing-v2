public class RayHit
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