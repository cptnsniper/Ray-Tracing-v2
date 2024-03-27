public class Vector2
{
    public double x;
    public double y;

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