package raytracer;

/**
 * Created by Jialin Gao on 3/31/2019.
 */

import org.joml.Vector4f;
import util.Material;

/**
 * This class Stores all the information that you will need to determine the closest object that was hit,
 * and information about that object to calculate shading.
 *
 * All the information contain
 * (a) the time ‘t’ on the ray where it intersects an object
 * (b) the 3D point of intersection in view coordinates
 * (c) the 3D normal of the object at that point in view coordinates
 * (d) The material properties.
 * (e) Texture coordinates and a Texture object, if applicable.
 */

public class HitRecord {
    private float t;
    private Vector4f intersectionP_view;
    private Vector4f intersectionP_obj;
    private Vector4f normal;
    private Material material;
    private boolean isHit;
    String textureName;
    float textureX;
    float textureY;

    // TODO: Texture

    /**
     * Constructors
     */

    public HitRecord() {
        t = 0;
        intersectionP_view = new Vector4f(0, 0, 0, 1);
        intersectionP_obj = new Vector4f(0, 0, 0, 1);
        normal = new Vector4f(0, 0, 0, 0);
        material = new Material();
        isHit = false;
    }

    public HitRecord(float t, Vector4f intersectionP_view, Vector4f intersectionP_obj, Vector4f normal, Material material) {
        this.t = t;
        this.intersectionP_view = new Vector4f(intersectionP_view);
        this.intersectionP_obj = new Vector4f(intersectionP_obj);
        this.normal = new Vector4f(normal);
        this.material = material;
        this.isHit = true;
    }

    /**
     * Setters and Getters
     */

    public float getT() {
        return t;
    }

    public void setT(float t) {
        this.t = t;
    }

    public Vector4f getNormal() {
        return normal;
    }

    public void setNormal(Vector4f normal) {
        this.normal = normal;
    }

    public Material getMaterial() {
        return material;
    }

    public void setMaterial(Material material) {
        this.material = material;
    }

    public boolean isHit() {
        return isHit;
    }

    public void setHit(boolean hit) {
        isHit = hit;
    }

    public Vector4f getIntersectionP_view() {
        return intersectionP_view;
    }

    public void setIntersectionP_view(Vector4f intersectionP_view) {
        this.intersectionP_view = intersectionP_view;
    }

    public Vector4f getIntersectionP_obj() {
        return intersectionP_obj;
    }

    public void setIntersectionP_obj(Vector4f intersectionP_obj) {
        this.intersectionP_obj = intersectionP_obj;
    }

    public String getTextureName() {
        return textureName;
    }

    public void setTextureName(String textureName) {
        this.textureName = textureName;
    }

    public float getTextureX() {
        return textureX;
    }

    public void setTextureX(float textureX) {
        this.textureX = textureX;
    }

    public float getTextureY() {
        return textureY;
    }

    public void setTextureY(float textureY) {
        this.textureY = textureY;
    }

}
