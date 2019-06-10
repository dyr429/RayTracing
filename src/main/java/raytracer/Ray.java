package raytracer;

import org.joml.Matrix4f;
import org.joml.Vector4f;

/**
 * Created by Jialin Gao on 3/31/2019.
 */

/**
 * This class contains a starting 3D point and a direction as a 3D vector.
 */

public class Ray {
    private Vector4f startPoint;
    private Vector4f direction;

    /**
     * Constructors
     */

    public Ray() {
        this.setStartPoint(new Vector4f(0, 0, 0, 1));
        this.setDirection(new Vector4f(0, 0, 0, 0));
    }

    public Ray(Vector4f startPoint, Vector4f direction) {
        this.setStartPoint(new Vector4f(startPoint));
        this.setDirection(new Vector4f(direction));
    }

    /**
     * Setters and Getters
     */

    public Vector4f getStartPoint() {
        return startPoint;
    }

    public void setStartPoint(Vector4f startPoint) {
        this.startPoint = new Vector4f(startPoint);
    }

    public Vector4f getDirection() {
        return direction;
    }

    public void setDirection(Vector4f direction) {
        this.direction = new Vector4f(direction);
    }

    /**
     * transform 3D ray to the world coordinate
     * @param modelview
     * @return a 3D ray in the
     */
    public Ray rayTransform(Matrix4f modelview) {
        Ray ray = new Ray(this.getStartPoint(), this.getDirection());
        ray.setStartPoint(ray.getStartPoint().mul(modelview));
        ray.setDirection(ray.getDirection().mul(modelview));
        return ray;
    }

    /**
     *
     * @param t the time on the ray where it intersects an object
     * @return the point of intersection
     */
    public Vector4f intersectObj(float t) {
        Vector4f point = new Vector4f(this.direction).mul(t);
        return new Vector4f(this.startPoint).add(point);
    }

}
