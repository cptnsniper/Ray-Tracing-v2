public class Vector3
{
    public double x;
    public double y;
    public double z;
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